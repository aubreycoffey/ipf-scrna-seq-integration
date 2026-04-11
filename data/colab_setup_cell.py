"""
================================================================================
IPF scRNA-seq Integration Project — Colab Setup & Data Loading
================================================================================
Paste this entire file as the first cell in each notebook when running in Colab.
For local / cloud VM use, set IN_COLAB = False and adjust PROJECT_DIR below.

What this does:
  1.  Detects environment (Colab vs local)
  2.  Mounts Google Drive (Colab only)
  3.  Installs required packages
  4.  Sets up project directory structure
  5.  Downloads GSE135893 + GSE136831 from NCBI GEO (skips files already present)
  6.  Loads each dataset into an AnnData object
  7.  Attaches donor / disease / cell-type metadata to each
  8.  Saves GSE135893_raw.h5ad and GSE136831_raw.h5ad to data/processed/
      for use by downstream notebooks

Approximate download sizes:
  GSE135893 : ~1.1 GB  (~5 min on a fast connection)
  GSE136831 : ~2.3 GB  (~10 min on a fast connection)
  Total disk : ~15 GB  (raw + processed .h5ad files)
================================================================================
"""

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
    # Edit this path if running locally
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
# The pandas version conflict warning is harmless — ignore it.

print("\nInstalling packages (this takes ~2 min on first run)...")
os.system(
    "pip install -q scanpy anndata harmonypy scvi-tools liana-py "
    "pyscenic sccoda gseapy GEOparse python-igraph leidenalg"
)
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

print(f"Scanpy  : {sc.__version__}")
print(f"AnnData : {ad.__version__}")

# ── 5. Shared download utility ────────────────────────────────────────────────

def download_with_progress(url: str, dest: str, label: str = "") -> None:
    """Download url → dest with a progress bar.
    Skips silently if dest already exists and is non-empty."""
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        print(f"  ✓ Already present : {os.path.basename(dest)} "
              f"({os.path.getsize(dest)/1e6:.1f} MB)")
        return

    print(f"  ↓ Downloading {os.path.basename(dest)} — {label}")

    def reporthook(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0:
            pct      = min(downloaded / total_size * 100, 100)
            mb_done  = downloaded / 1e6
            mb_total = total_size / 1e6
            filled   = int(30 * pct / 100)
            bar      = "█" * filled + "░" * (30 - filled)
            print(f"\r    [{bar}] {pct:5.1f}%  {mb_done:.0f}/{mb_total:.0f} MB",
                  end="", flush=True)

    urllib.request.urlretrieve(url, dest, reporthook)
    print(f"\r    ✓ Done — {os.path.getsize(dest)/1e6:.1f} MB")


def check_downloads(dest_dir: str, filenames: list) -> None:
    """Verify all files downloaded successfully. Raises if any are missing."""
    print("\nFile check:")
    all_ok = True
    for fname in filenames:
        path = f"{dest_dir}/{fname}"
        if os.path.exists(path) and os.path.getsize(path) > 1000:
            print(f"  ✓ {fname} ({os.path.getsize(path)/1e6:.1f} MB)")
        else:
            print(f"  ✗ MISSING or EMPTY : {fname}")
            all_ok = False
    if not all_ok:
        raise RuntimeError(
            "One or more files failed to download. "
            "Re-run this cell — files already present will be skipped."
        )


# ═════════════════════════════════════════════════════════════════════════════
# DATASET 1 — GSE135893  (Habermann et al. 2020)
# 114,396 cells | 20 PF + 10 controls | Vanderbilt / TGen cohort
# ═════════════════════════════════════════════════════════════════════════════

print("\n" + "═"*72)
print("DATASET 1 — GSE135893  (Habermann et al. 2020)")
print("═"*72)

GSE1_DIR  = f"{RAW_DIR}/GSE135893"
GSE1_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135893/suppl"
os.makedirs(GSE1_DIR, exist_ok=True)

GSE1_FILES = {
    "GSE135893_matrix.mtx.gz"      : "count matrix (~1 GB)",
    "GSE135893_barcodes.tsv.gz"     : "cell barcodes (~1 MB)",
    "GSE135893_genes.tsv.gz"        : "gene names (~0.1 MB)",
    "GSE135893_IPF_metadata.csv.gz" : "donor metadata (~3 MB)",
}

print("\n── Downloading ────────────────────────────────────────────────────────")
for fname, label in GSE1_FILES.items():
    download_with_progress(f"{GSE1_BASE}/{fname}", f"{GSE1_DIR}/{fname}", label)

check_downloads(GSE1_DIR, list(GSE1_FILES.keys()))

# ── Load GSE135893 ────────────────────────────────────────────────────────────

print("\n── Loading ────────────────────────────────────────────────────────────")

print("  Reading barcodes...")
with gzip.open(f"{GSE1_DIR}/GSE135893_barcodes.tsv.gz", "rt") as f:
    barcodes1 = [line.strip() for line in f]
print(f"  {len(barcodes1):,} cells")

print("  Reading genes...")
with gzip.open(f"{GSE1_DIR}/GSE135893_genes.tsv.gz", "rt") as f:
    genes_raw = [line.strip().split("\t") for line in f]
if len(genes_raw[0]) >= 2:
    gene_ids1     = [g[0] for g in genes_raw]
    gene_symbols1 = [g[1] for g in genes_raw]
else:
    gene_ids1     = [g[0] for g in genes_raw]
    gene_symbols1 = gene_ids1
print(f"  {len(gene_ids1):,} genes  |  example: {gene_symbols1[:3]}")

print("  Reading count matrix (~1–2 min)...")
with gzip.open(f"{GSE1_DIR}/GSE135893_matrix.mtx.gz", "rb") as f:
    X1 = scipy.io.mmread(f).T.tocsr()
print(f"  Matrix shape: {X1.shape}  (cells × genes)")

print("  Assembling AnnData...")
adata1 = ad.AnnData(
    X   = X1,
    obs = pd.DataFrame(index=barcodes1),
    var = pd.DataFrame(
        {"gene_ids": gene_ids1, "gene_symbols": gene_symbols1},
        index=gene_symbols1,
    ),
)
adata1.var_names_make_unique()

# Attach metadata
meta1 = pd.read_csv(f"{GSE1_DIR}/GSE135893_IPF_metadata.csv.gz", index_col=0)
print(f"\n  Metadata shape   : {meta1.shape}")
print(f"  Metadata columns : {list(meta1.columns)}")

shared1 = adata1.obs.index.intersection(meta1.index)
if len(shared1) == 0:
    adata1.obs.index = adata1.obs.index.str.replace(r"-\d+$", "", regex=True)
    shared1 = adata1.obs.index.intersection(meta1.index)
print(f"  Barcodes matched : {len(shared1):,} / {adata1.n_obs:,}")

if len(shared1) > 0:
    adata1.obs = adata1.obs.join(meta1, how="left")
else:
    print("  WARNING: barcodes could not be matched to metadata.")
    print(f"    adata  (first 3): {list(adata1.obs.index[:3])}")
    print(f"    meta   (first 3): {list(meta1.index[:3])}")

adata1.obs["dataset"] = "GSE135893"

print(f"\n  adata1 : {adata1}")

# Save
SAVE1 = f"{PROCESSED_DIR}/GSE135893_raw.h5ad"
print(f"\n  Saving → {SAVE1}")
adata1.write_h5ad(SAVE1, compression="gzip")
print(f"  ✓ Saved ({os.path.getsize(SAVE1)/1e6:.0f} MB)")


# ═════════════════════════════════════════════════════════════════════════════
# DATASET 2 — GSE136831  (Adams et al. 2020)
# 312,928 cells | 32 IPF + 18 COPD + 28 controls | Yale / Harvard cohort
# ═════════════════════════════════════════════════════════════════════════════

print("\n" + "═"*72)
print("DATASET 2 — GSE136831  (Adams et al. 2020)")
print("═"*72)

GSE2_DIR  = f"{RAW_DIR}/GSE136831"
GSE2_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136831/suppl"
os.makedirs(GSE2_DIR, exist_ok=True)

# Note: GSE136831_RAW.tar is 16 GB and contains per-sample files we don't need.
# We use the pre-assembled combined matrix files instead.
GSE2_FILES = {
    "GSE136831_RawCounts_Sparse.mtx.gz"                    : "count matrix (~2 GB)",
    "GSE136831_AllCells.cellBarcodes.txt.gz"                : "cell barcodes (~1.3 MB)",
    "GSE136831_AllCells.GeneIDs.txt.gz"                     : "gene names (~316 KB)",
    "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz": "metadata incl. disease + cell type (~4 MB)",
}

print("\n── Downloading ────────────────────────────────────────────────────────")
for fname, label in GSE2_FILES.items():
    download_with_progress(f"{GSE2_BASE}/{fname}", f"{GSE2_DIR}/{fname}", label)

check_downloads(GSE2_DIR, list(GSE2_FILES.keys()))

# ── Load GSE136831 ────────────────────────────────────────────────────────────

print("\n── Loading ────────────────────────────────────────────────────────────")

print("  Reading barcodes...")
with gzip.open(f"{GSE2_DIR}/GSE136831_AllCells.cellBarcodes.txt.gz", "rt") as f:
    lines2 = [line.strip().replace('"', '') for line in f if line.strip()]
# Skip header row if present
if lines2 and any(k in lines2[0].upper() for k in ["BARCODE", "CELL", "ID"]):
    lines2 = lines2[1:]
barcodes2 = lines2
print(f"  {len(barcodes2):,} cells")

print("  Reading genes...")
with gzip.open(f"{GSE2_DIR}/GSE136831_AllCells.GeneIDs.txt.gz", "rt") as f:
    genes_raw2 = [line.strip().replace('"', '').split("\t") for line in f]

# GSE136831 has a header row — skip it if first entry looks like a column name
first = genes_raw2[0][0].upper()
if any(k in first for k in ["ID", "GENE", "HGNC", "ENSEMBL", "SYMBOL"]):
    genes_raw2 = genes_raw2[1:]

# Determine column layout: [ensembl_id, gene_symbol] or [gene_symbol] only
if len(genes_raw2[0]) >= 2:
    gene_ids2     = [g[0] for g in genes_raw2]
    gene_symbols2 = [g[1] for g in genes_raw2]
else:
    gene_ids2     = [g[0] for g in genes_raw2]
    gene_symbols2 = gene_ids2

print(f"  {len(gene_ids2):,} genes  |  example: {gene_symbols2[:3]}")

print("  Reading count matrix (~3–5 min for 2 GB)...")
with gzip.open(f"{GSE2_DIR}/GSE136831_RawCounts_Sparse.mtx.gz", "rb") as f:
    X2 = scipy.io.mmread(f).T.tocsr()
print(f"  Matrix shape: {X2.shape}  (cells × genes)")

print("  Assembling AnnData...")
adata2 = ad.AnnData(
    X   = X2,
    obs = pd.DataFrame(index=barcodes2),
    var = pd.DataFrame(
        {"gene_ids": gene_ids2, "gene_symbols": gene_symbols2},
        index=gene_symbols2,
    ),
)
adata2.var_names_make_unique()

# Attach metadata
# GSE136831 metadata is rich — includes Sample_Name, Disease, CellType_Category,
# CellType_Anno, donor info. We read it and join on barcode.
meta2 = pd.read_csv(
    f"{GSE2_DIR}/GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz",
    sep="\t",
    index_col=0,
)
print(f"\n  Metadata shape   : {meta2.shape}")
print(f"  Metadata columns : {list(meta2.columns)}")
print(f"\n  Disease groups:")
if "Disease_Identity" in meta2.columns:
    print(meta2["Disease_Identity"].value_counts().to_string())
elif "Status" in meta2.columns:
    print(meta2["Status"].value_counts().to_string())
else:
    # Print first column that looks like disease status
    for col in meta2.columns:
        if any(k in col.lower() for k in ["disease", "status", "diagnosis"]):
            print(meta2[col].value_counts().to_string())
            break

shared2 = adata2.obs.index.intersection(meta2.index)
if len(shared2) == 0:
    adata2.obs.index = adata2.obs.index.str.replace(r"-\d+$", "", regex=True)
    shared2 = adata2.obs.index.intersection(meta2.index)
print(f"\n  Barcodes matched : {len(shared2):,} / {adata2.n_obs:,}")

if len(shared2) > 0:
    adata2.obs = adata2.obs.join(meta2, how="left")
else:
    print("  WARNING: barcodes could not be matched to metadata.")
    print(f"    adata  (first 3): {list(adata2.obs.index[:3])}")
    print(f"    meta   (first 3): {list(meta2.index[:3])}")

adata2.obs["dataset"] = "GSE136831"

print(f"\n  adata2 : {adata2}")

# Save
SAVE2 = f"{PROCESSED_DIR}/GSE136831_raw.h5ad"
print(f"\n  Saving → {SAVE2}")
adata2.write_h5ad(SAVE2, compression="gzip")
print(f"  ✓ Saved ({os.path.getsize(SAVE2)/1e6:.0f} MB)")


# ── Sanity check — both datasets ─────────────────────────────────────────────

print("\n" + "═"*72)
print("SANITY CHECK")
print("═"*72)

markers = ["SFTPC", "SFTPB", "AGER", "KRT17", "ACTA2",
           "DCN", "SPP1", "FOXJ1", "SCGB1A1", "CD14", "COL1A1", "CA4"]

for label, adata in [("GSE135893", adata1), ("GSE136831", adata2)]:
    found   = [g for g in markers if g in adata.var_names]
    missing = [g for g in markers if g not in adata.var_names]
    print(f"\n  {label}")
    print(f"    Cells   : {adata.n_obs:,}")
    print(f"    Genes   : {adata.n_vars:,}")
    print(f"    Markers found   : {found}")
    if missing:
        print(f"    Markers missing : {missing}  (a few missing is normal)")


# ── Summary ───────────────────────────────────────────────────────────────────

print(f"""
{"═"*72}
  Setup complete. Both datasets downloaded, loaded, and saved.

  adata1  (GSE135893) : {adata1.n_obs:,} cells × {adata1.n_vars:,} genes
  adata2  (GSE136831) : {adata2.n_obs:,} cells × {adata2.n_vars:,} genes

  Saved to:
    data/processed/GSE135893_raw.h5ad
    data/processed/GSE136831_raw.h5ad

  In subsequent notebooks, reload with:
    adata1 = sc.read_h5ad(f"{{PROCESSED_DIR}}/GSE135893_raw.h5ad")
    adata2 = sc.read_h5ad(f"{{PROCESSED_DIR}}/GSE136831_raw.h5ad")

  Next step: 01_qc_preprocessing.ipynb
{"═"*72}
""")

# ── Subsample flag (for fast development runs) ────────────────────────────────
# Set SUBSAMPLE = True while building / testing notebooks.
# Set back to False before running the full analysis.

SUBSAMPLE   = False
SUBSAMPLE_N = 20_000

if SUBSAMPLE:
    sc.pp.subsample(adata1, n_obs=SUBSAMPLE_N, random_state=42)
    sc.pp.subsample(adata2, n_obs=SUBSAMPLE_N, random_state=42)
    print(f"  Subsampled both datasets to {SUBSAMPLE_N:,} cells for development.")
