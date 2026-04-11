"""
================================================================================
IPF scRNA-seq Integration Project — Colab Setup & Data Loading
================================================================================
Run this file to set up the data in Google Drive. This only needs to be done 
once. For subsequent notebooks, just run the header at the top of this file. 

For local / cloud VM use, set IN_COLAB = False and adjust PROJECT_DIR below.

What this does:
  1. Detects environment (Colab vs local)
  2. Mounts Google Drive (Colab only)
  3. Installs required packages
  4. Sets up project directory structure
  5. Downloads GSE135893 from NCBI GEO (skips files already present)
  6. Loads GSE135893 into an AnnData object
  7. Attaches donor/disease metadata
  8. Saves a processed .h5ad to data/processed/ for downstream notebooks
================================================================================
"""



# Standard header for notebooks 01–07
import sys, os
IN_COLAB = "google.colab" in sys.modules
if IN_COLAB:
    from google.colab import drive
    drive.mount("/content/drive")
    PROJECT_DIR = "/content/drive/MyDrive/ipf-scrna-seq-integration"
else:
    PROJECT_DIR = os.path.dirname(os.path.abspath("__file__"))

DATA_DIR      = f"{PROJECT_DIR}/data"
PROCESSED_DIR = f"{DATA_DIR}/processed"
RESULTS_DIR   = f"{PROJECT_DIR}/results/figures"

import scanpy as sc
adata = sc.read_h5ad(f"{PROCESSED_DIR}/GSE135893_raw.h5ad")
print(adata)







# ── 0. Imports available before any installs ─────────────────────────────────

import sys
import os

# ── 1. Environment detection ──────────────────────────────────────────────────

IN_COLAB = "google.colab" in sys.modules

print(f"Environment: {'Google Colab' if IN_COLAB else 'Local / Cloud VM'}")

# ── 2. Mount Google Drive (Colab only) ───────────────────────────────────────

if IN_COLAB:
    from google.colab import drive
    drive.mount("/content/drive")
    PROJECT_DIR = "/content/drive/MyDrive/ipf-scrna-seq-integration"
else:
    # ── Edit this path if running locally ────────────────────────────────────
    PROJECT_DIR = os.path.dirname(os.path.abspath("__file__"))

DATA_DIR      = f"{PROJECT_DIR}/data"
RAW_DIR       = f"{DATA_DIR}/raw"
PROCESSED_DIR = f"{DATA_DIR}/processed"
RESULTS_DIR   = f"{PROJECT_DIR}/results/figures"

for d in [RAW_DIR, PROCESSED_DIR, RESULTS_DIR]:
    os.makedirs(d, exist_ok=True)

print(f"Project dir : {PROJECT_DIR}")
print(f"Data dir    : {DATA_DIR}")

# ── 3. Install packages ───────────────────────────────────────────────────────
# Suppress the pandas version conflict warning — it is harmless.

print("\nInstalling packages (this takes ~2 min on first run)...")

os.system("pip install -q scanpy anndata harmonypy scvi-tools liana-py pyscenic sccoda gseapy GEOparse python-igraph leidenalg")

print("Packages installed.")

# ── 4. Imports (post-install) ─────────────────────────────────────────────────

import urllib.request
import gzip
import scipy.io
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc

sc.settings.verbosity = 2
sc.settings.figdir    = RESULTS_DIR
sc.settings.set_figure_params(dpi=100, facecolor="white")

print(f"Scanpy version : {sc.__version__}")
print(f"AnnData version: {ad.__version__}")

# ── 5. Download GSE135893 from NCBI GEO ──────────────────────────────────────

GSE1_DIR  = f"{RAW_DIR}/GSE135893"
GSE1_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl"

os.makedirs(GSE1_DIR, exist_ok=True)

GSE1_FILES = {
    "GSE135893_matrix.mtx.gz"      : "count matrix (~1 GB, takes ~5 min)",
    "GSE135893_barcodes.tsv.gz"     : "cell barcodes (~1 MB)",
    "GSE135893_genes.tsv.gz"        : "gene names (~0.1 MB)",
    "GSE135893_IPF_metadata.csv.gz" : "donor metadata (~3 MB)",
}

def download_with_progress(url: str, dest: str, label: str = "") -> None:
    """Download url to dest with progress bar. Skips if file already exists."""
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"  ✓ Already present: {os.path.basename(dest)} "
              f"({os.path.getsize(dest)/1e6:.1f} MB)")
        return

    print(f"  ↓ Downloading {os.path.basename(dest)} — {label}")

    def reporthook(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0:
            pct      = min(downloaded / total_size * 100, 100)
            mb_done  = downloaded / 1e6
            mb_total = total_size / 1e6
            bar_len  = 30
            filled   = int(bar_len * pct / 100)
            bar      = "█" * filled + "░" * (bar_len - filled)
            print(f"\r    [{bar}] {pct:5.1f}%  {mb_done:.0f}/{mb_total:.0f} MB",
                  end="", flush=True)

    urllib.request.urlretrieve(url, dest, reporthook)
    print(f"\r    ✓ Done — {os.path.getsize(dest)/1e6:.1f} MB")

print("\n── Downloading GSE135893 ──────────────────────────────────────────────")
for fname, label in GSE1_FILES.items():
    download_with_progress(f"{GSE1_BASE}/{fname}", f"{GSE1_DIR}/{fname}", label)

# Verify all files are present and non-empty
print("\nFile check:")
all_ok = True
for fname in GSE1_FILES:
    path = f"{GSE1_DIR}/{fname}"
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        print(f"  ✓ {fname} ({os.path.getsize(path)/1e6:.1f} MB)")
    else:
        print(f"  ✗ MISSING or EMPTY: {fname}")
        all_ok = False

if not all_ok:
    raise RuntimeError(
        "One or more files failed to download. "
        "Re-run this cell — download_with_progress skips files already present."
    )

# ── 6. Load count matrix into AnnData ─────────────────────────────────────────

print("\n── Loading count matrix ───────────────────────────────────────────────")

# barcodes (cell IDs)
print("  Reading barcodes...")
with gzip.open(f"{GSE1_DIR}/GSE135893_barcodes.tsv.gz", "rt") as f:
    barcodes = [line.strip() for line in f]
print(f"  {len(barcodes):,} cells")

# genes
print("  Reading genes...")
with gzip.open(f"{GSE1_DIR}/GSE135893_genes.tsv.gz", "rt") as f:
    genes_raw = [line.strip().split("\t") for line in f]

# GEO gene files can be [ensembl_id] or [ensembl_id, gene_symbol]
if len(genes_raw[0]) >= 2:
    gene_ids     = [g[0] for g in genes_raw]
    gene_symbols = [g[1] for g in genes_raw]
else:
    gene_ids     = [g[0] for g in genes_raw]
    gene_symbols = gene_ids

print(f"  {len(gene_ids):,} genes  |  example: {gene_symbols[:3]}")

# count matrix
print("  Reading count matrix (this takes ~1–2 min)...")
with gzip.open(f"{GSE1_DIR}/GSE135893_matrix.mtx.gz", "rb") as f:
    X = scipy.io.mmread(f).T.tocsr()   # transpose → rows=cells, cols=genes

print(f"  Matrix shape: {X.shape}  (cells × genes)")

# assemble AnnData
print("  Assembling AnnData...")
adata = ad.AnnData(
    X   = X,
    obs = pd.DataFrame(index=barcodes),
    var = pd.DataFrame(
        {"gene_ids": gene_ids, "gene_symbols": gene_symbols},
        index=gene_symbols,
    ),
)
adata.var_names_make_unique()

print(f"\n  {adata}")

# ── 7. Attach donor / disease metadata ────────────────────────────────────────

print("\n── Attaching metadata ─────────────────────────────────────────────────")

meta = pd.read_csv(f"{GSE1_DIR}/GSE135893_IPF_metadata.csv.gz", index_col=0)
print(f"  Metadata shape : {meta.shape}")
print(f"  Columns        : {list(meta.columns)}")
print(f"\n  Preview:\n{meta.head()}")

# Join by barcode — try as-is first, then strip trailing -1 suffix
shared = adata.obs.index.intersection(meta.index)
if len(shared) == 0:
    adata.obs.index = adata.obs.index.str.replace(r"-\d+$", "", regex=True)
    shared = adata.obs.index.intersection(meta.index)

print(f"\n  Barcodes matched: {len(shared):,} / {adata.n_obs:,}")

if len(shared) > 0:
    adata.obs = adata.obs.join(meta, how="left")
    for col in meta.columns:
        print(f"    {col}: {adata.obs[col].nunique()} unique values — "
              f"{list(adata.obs[col].value_counts().head(3).index)}")
else:
    print("  WARNING: barcodes could not be matched to metadata.")
    print(f"    adata barcodes (first 3) : {list(adata.obs.index[:3])}")
    print(f"    metadata index  (first 3): {list(meta.index[:3])}")

adata.obs["dataset"] = "GSE135893"

# ── 8. Sanity check ───────────────────────────────────────────────────────────

print("\n── Sanity check ───────────────────────────────────────────────────────")

markers = ["SFTPC", "SFTPB", "AGER", "KRT17", "ACTA2", "DCN",
           "SPP1", "FOXJ1", "SCGB1A1", "CD14", "COL1A1", "CA4"]
found   = [g for g in markers if g in adata.var_names]
missing = [g for g in markers if g not in adata.var_names]

print(f"  Marker genes found   : {found}")
if missing:
    print(f"  Marker genes missing : {missing}  (a few missing is normal)")

print(f"\n  Cells : {adata.n_obs:,}")
print(f"  Genes : {adata.n_vars:,}")

# ── 9. Save raw AnnData ───────────────────────────────────────────────────────

SAVE_PATH = f"{PROCESSED_DIR}/GSE135893_raw.h5ad"
print(f"\n── Saving → {SAVE_PATH}")
adata.write_h5ad(SAVE_PATH, compression="gzip")
print(f"  ✓ Saved ({os.path.getsize(SAVE_PATH)/1e6:.0f} MB)")

# ── 10. Summary ───────────────────────────────────────────────────────────────

print(f"""
════════════════════════════════════════════════════════════════════════════════
  Setup complete.

  adata                    : {adata.n_obs:,} cells × {adata.n_vars:,} genes
  Saved to                 : data/processed/GSE135893_raw.h5ad

  In subsequent notebooks, reload with:
      adata = sc.read_h5ad(f"{{PROCESSED_DIR}}/GSE135893_raw.h5ad")

  Next step: 01_qc_preprocessing.ipynb
════════════════════════════════════════════════════════════════════════════════
""")

# ── 11. Subsample flag ────────────────────────────────────────────────────────
# Set SUBSAMPLE = True to work on a 20k-cell slice during development.
# Set back to False before running the full analysis.

SUBSAMPLE   = False
SUBSAMPLE_N = 20_000

if SUBSAMPLE:
    sc.pp.subsample(adata, n_obs=SUBSAMPLE_N, random_state=42)
    print(f"  Subsampled to {adata.n_obs:,} cells for development.")
