# ipf-scrna-seq-integration

**Cross-cohort integrative single-cell transcriptomic analysis of idiopathic pulmonary fibrosis: trajectory inference, cell-cell communication, and transcription factor regulon analysis across all major lung compartments.**

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
│   ├── 02_integration_clustering.ipynb
│   ├── 03_differential_analysis.ipynb
│   ├── 04_trajectory_inference.ipynb
│   ├── 05_cellchat_communication.ipynb
│   ├── 06_pyscenic_regulons.ipynb
│   └── 07_validation.ipynb
│
├── src/
│   ├── preprocessing.py
│   ├── integration.py
│   ├── annotation.py
│   ├── trajectory.py
│   ├── communication.py
│   ├── regulons.py
│   └── plotting.py
│
└── results/
    ├── figures/
    └── tables/
```

---

## Reproducing the Analysis

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/ipf-scrna-seq-integration.git
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
02 → Cross-cohort integration, UMAP, clustering, annotation
03 → Differential abundance and expression (IPF vs. control vs. COPD)
04 → Trajectory inference (epithelial, mesenchymal, immune)
05 → Cell-cell communication (LIANA)
06 → Transcription factor regulon analysis (pySCENIC)
07 → Cross-cohort validation and COPD specificity
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

## Key Findings

*(To be updated upon preprint submission)*

---

## Citation

If you use this code or analysis, please cite:

```bibtex
@article{YOUR_NAME_ipf_2025,
  title   = {Cross-cohort integrative single-cell transcriptomic analysis reveals
             coordinated cell state transitions, intercellular signaling rewiring,
             and transcription factor programs in idiopathic pulmonary fibrosis},
  author  = {YOUR NAME},
  journal = {arXiv},
  year    = {2025},
  url     = {https://arxiv.org/abs/XXXX.XXXXX}
}
```

Please also cite the original data sources:

- Habermann et al. (2020) *Science Advances* doi:10.1126/sciadv.aba1972
- Adams et al. (2020) *Science Advances* doi:10.1126/sciadv.aba1983

---

## License

MIT License. See [`LICENSE`](LICENSE) for details.

---

## Contact

Questions or issues? Open a GitHub issue or reach out at aubreymariecoffey@gmail.com.
