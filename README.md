# ipf-scrna-seq-integration

**Cross-cohort integrative single-cell transcriptomic analysis reveals coordinated cell state transitions, intercellular signaling rewiring, and transcription factor programs in idiopathic pulmonary fibrosis.**

---

## Overview

Idiopathic pulmonary fibrosis (IPF) is a fatal progressive lung disease with a median survival of 3–5 years. Two landmark single-cell RNA sequencing (scRNA-seq) studies — Habermann et al. 2020 and Adams et al. 2020 — characterized the major disease-associated cell populations in independent cohorts. However, no study has integrated both cohorts, performed systematic cell-cell communication analysis across all major cell types, or extended transcription factor regulon analysis beyond the epithelial compartment.

This project addresses those gaps. We integrate two large publicly available IPF scRNA-seq datasets (>400,000 cells across 52 IPF, 46 control, and 18 COPD donors) and apply an end-to-end Python/Scanpy pipeline encompassing:

- **Cross-cohort batch correction** (Harmony + scVI) and unified cell type annotation
- **Multi-compartment trajectory inference** (PAGA + diffusion pseudotime) across epithelial, mesenchymal, and immune lineages — including a direct cross-cohort test of Adams et al.'s contested finding of limited fibroblast-to-myofibroblast connectivity
- **Systematic cell-cell communication analysis** using LIANA (aggregating CellChat, CellPhoneDB, and NATMI), comparing IPF vs. control vs. COPD — the primary novel contribution of this work
- **Transcription factor regulon analysis** (pySCENIC) extended to mesenchymal, immune, and endothelial compartments, beyond the epithelial-only analysis of prior work

---

## Datasets

| Dataset | GEO Accession | Cells | Donors | Reference |
|---|---|---|---|---|
| Habermann et al. 2020 | [GSE135893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893) | ~114,000 | 20 PF + 10 control | *Science Advances* 2020 |
| Adams et al. 2020 | [GSE136831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831) | ~312,000 | 32 IPF + 18 COPD + 28 control | *Science Advances* 2020 |

Raw data is not stored in this repository. See [`data/README.md`](data/README.md) for download instructions or run:

```bash
bash data/download_data.sh
```

---

## Repository Structure

```
ipf-scrna-seq-integration/
│
├── README.md
├── environment.yml             
├── requirements.txt
│
├── data/
│   ├── README.md
│   │── colab_setup_cell.py               
│   └── download_data.sh        
│
├── notebooks/
│   ├── 01_qc_preprocessing.ipynb

```

---

## Reproducing the Analysis

### 1. Clone the repository

```bash
git clone https://github.com/aubreycoffey/ipf-scrna-seq-integration.git
cd ipf-scrna-seq-integration
```

### 2. Set up the environment

```bash
conda env create -f environment.yml
conda activate ipf-scrna
```

### 3. Download data

```bash
bash data/download_data.sh
```

### 4. Run notebooks in order

Notebooks are numbered and should be run sequentially. Each notebook saves intermediate `.h5ad` files to `data/processed/` for use by downstream notebooks.

```
01 → QC and preprocessing (per dataset)
```

---

## Key Methods

| Step | Tool | Notes |
|---|---|---|
| QC & normalization | Scanpy | Per-dataset, then merged |
| Batch correction | Harmony, scVI | Corrects for donor + dataset effects |
| Clustering | Leiden | Multiple resolutions evaluated |
| Trajectory inference | PAGA + diffusion pseudotime | Scanpy |
| Cell-cell communication | LIANA (liana-py) | Aggregates CellChat, CellPhoneDB, NATMI |
| Regulon analysis | pySCENIC | GRNBoost2 → cisTarget → AUCell |
| Compositional testing | scCODA | Bayesian; IPF vs. control vs. COPD |
| GSEA | GSEApy | Hallmark + KEGG gene sets |


---

## Contact

Questions or issues? Open a GitHub issue or reach out at aubreymariecoffey@gmail.com.
