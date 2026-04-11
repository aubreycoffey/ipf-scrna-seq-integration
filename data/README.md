# Data

Raw data is not stored in this repository. Both datasets are publicly available
on NCBI GEO and licensed under CC BY 4.0.

---

## Datasets

### GSE135893 — Habermann et al. 2020
- **Paper:** Single-cell RNA sequencing reveals profibrotic roles of distinct epithelial
  and mesenchymal lineages in pulmonary fibrosis. *Science Advances* 2020.
- **GEO page:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893
- **Cells:** ~114,000
- **Donors:** 20 PF (12 IPF + 8 other ILD) + 10 non-fibrotic controls
- **Format:** 10x Genomics count matrices (barcodes, features, matrix)
- **Size on disk:** ~3 GB compressed

### GSE136831 — Adams et al. 2020
- **Paper:** Single-cell RNA-seq reveals ectopic and aberrant lung-resident cell
  populations in idiopathic pulmonary fibrosis. *Science Advances* 2020.
- **GEO page:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831
- **Cells:** ~312,000
- **Donors:** 32 IPF + 18 COPD + 28 controls
- **Format:** 10x Genomics count matrices per sample
- **Size on disk:** ~8 GB compressed

**Total disk requirement (raw):** ~11 GB compressed, ~25 GB uncompressed  
**Peak RAM during analysis:** ~40–60 GB (full integrated dataset)  
**Recommended minimum RAM:** 32 GB (you can subset to one dataset for development)

---

## Download Options

### Option 1: Automated download script (recommended)

```bash
bash data/download_data.sh
```

This script downloads both datasets using the GEO FTP server. It will create the
following directory structure:

```
data/
├── raw/
│   ├── GSE135893/
│   │   └── [sample directories with matrix files]
│   └── GSE136831/
│       └── [sample directories with matrix files]
├── processed/        # created by notebooks; .h5ad files saved here
└── README.md
```

### Option 2: Manual download from GEO

1. Go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893
2. Scroll to the bottom and click "Download family" → "SOFT formatted family file(s)"
   OR use the supplementary files section to download count matrices directly
3. Repeat for GSE136831

### Option 3: Download via Python (inside a notebook)

```python
import GEOparse

# Download GSE135893
gse1 = GEOparse.get_GEO(geo="GSE135893", destdir="data/raw/")

# Download GSE136831
gse2 = GEOparse.get_GEO(geo="GSE136831", destdir="data/raw/")
```

### Option 4: Google Colab / cloud environment

If you are running this project in Colab or a cloud VM, use the provided
Colab setup cell at the top of each notebook, which mounts Google Drive and
downloads only the files needed for that notebook.

---

## Storage Tips for Resource-Constrained Environments

- Process one dataset at a time and save the resulting `.h5ad` to persistent storage
  (Google Drive, S3, etc.) before loading the second
- The processed `.h5ad` files (post-QC, normalized, with embeddings) are much smaller
  than the raw count matrices — typically 1–3 GB each
- For development and testing, use the `--n_cells` flag in `src/preprocessing.py`
  to subsample to 10,000–20,000 cells before running the full pipeline
- pySCENIC (notebook 06) is the most memory-intensive step — run it on the full
  dataset only after confirming the pipeline works on a subsample

---

## Processed Data

Intermediate `.h5ad` files produced by each notebook are saved to `data/processed/`
and are used as inputs by downstream notebooks:

| File | Produced by | Used by |
|------|-------------|---------|
| `gse135893_qc.h5ad` | 01 | 02 |
| `gse136831_qc.h5ad` | 01 | 02 |
| `integrated.h5ad` | 02 | 03, 04, 05, 06, 07 |
| `integrated_with_trajectories.h5ad` | 04 | 06, 07 |
| `liana_results.h5ad` | 05 | 07 |
| `pyscenic_results.h5ad` | 06 | 07 |

These processed files are not stored in the repository but can be regenerated
by running the notebooks in order.
