# Aging signatures in ICE-aged brain organoids vs. the aging human brain

**link for GitHub Pages https://vkorobeynyk.github.io/aging_human_organoid_comparisson/**

Analysis code for comparing transcriptomic aging signatures between **ICE-treated ("aged") human brain organoids** and **natural human brain aging** using single-cell/single-nucleus RNA-seq.

ICE (Inducible Changes to the Epigenome; [Yang et al., Cell 2023](https://www.cell.com/cell/fulltext/S0092-8674(22)01570-7)) is used to induce an aged-like state in organoids. The central question is whether the gene-expression changes induced by ICE in organoids recapitulate those that occur with chronological age in human neurons. To address this, each dataset is analysed separately (QC, differential expression, GSEA, glmnet models), and the results are then compared across datasets — including training models on one dataset and testing them on the other.

## Contact

Vladyslav Korobeynyk, HIFO / DMLS, University of Zurich (Jessberger lab / Mark D. Robinson lab)

## Generative AI statement
Generative AI was used throughout this entire benchmark to make code nicer to read and more efficient. The entire logic of the benchmark was created by myself and I assume responsability of the content within this repo.

## Repository structure

```
aging_human_organoid_comparisson/
├── helper_functions.R                    # shared plotting / DEA / GSEA / modelling helpers
├── general_qc_sc_pipeline_R_function.R   # shared single-cell QC & processing functions
├── environment.yml                       # conda environment
├── LICENSE                               # GPL-3.0 licence
├── move_reports_to_docs.sh               # moves rendered .html reports into docs/
├── docs/                                 # rendered .html reports (GitHub Pages)
├── Organoid_analysis/
│   ├── 01_organoid_QC_processing_integration.qmd
│   ├── 02_organoid_Analysis_metacell_Level.qmd
│   ├── 03_organoid_Analysis_singleCell_Level.qmd
│   ├── 04_organoid_ML_predictions_organoids.qmd
│   ├── sbatch_analysis_organoids
│   ├── logs/                             # SLURM logs (contents not tracked)
│   └── organoid_analysis_output/         # DEA, GSEA and glmnet result tables
├── Human_analysis/
│   ├── 01_human_Processing.qmd
│   ├── 02_human_Analysis_metacell_Level.qmd
│   ├── 02b_human_Analysis_singleCell_Level.qmd
│   ├── 03_human_ML_predictions_human.qmd
│   ├── 04_human_ML_predictions_human_testOrganoids.qmd
│   ├── sbatch_analysis_human
│   ├── logs/                             # SLURM logs (contents not tracked)
│   └── human_analysis_output/            # DEA, GSEA and glmnet result tables
├── compare_organoid_human/
│   ├── 01_Comparison_human_organoid.qmd
│   └── sbatch_comparison_analysis
└── publication_plots/                    # figures exported for the manuscript
```

Each `.qmd` starts with a **"What this script does"** callout summarising its inputs, methods and outputs. Rendered `.html` reports are stored in `docs/` (same sub-folder structure) and listed in `docs/index.html`.

## Data

| Dataset | Description | Source |
|---|---|---|
| Organoids | 10x scRNA-seq of brain organoids, wild-type vs. ICE, across several dataset groups/timepoints (COp1 10w, COp1 4m, COp2 4dpi, COp2 4wpi) | In-house (not public yet) |
| Human | snRNA-seq of human brain, restricted to inhibitory and excitatory neurons of patients < 80 years with ≥ 500 cells per cell type | Raw: [EGA EGAS00001006345](https://ega-archive.org/studies/EGAS00001006345) · Processed: [Zenodo 13343729](https://zenodo.org/records/13343729) · Paper: [Nature Genetics (2024)](https://www.nature.com/articles/s41588-024-02050-9) |
| iNeurons | Bulk RNA-seq of induced neurons directly reprogrammed from fibroblasts of donors of different ages (22 runs, ERR668328–ERR668411), used as an additional test set for the human age model. FASTQs were processed to a count table with [ARMOR](https://github.com/csoneson/ARMOR) ([Orjuela et al., 2019](https://doi.org/10.1534/g3.119.400185)) and stored as `edgeR_dge_induced_Neurons.rds` | Paper: [Mertens et al., Cell Stem Cell (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5929130/) |

Raw data, `.rds` objects and job logs are **not tracked** in git (see `.gitignore`). They are expected inside the project folder (`organoid_data/`, `human_data/`, etc.).

### Abbreviations

| Term | Meaning |
|---|---|
| COp1 / COp2 | Cortical organoids generated with protocol 1 / protocol 2 |
| 10w / 4m | 10-week-old / 4-month-old organoids (COp1) |
| 4dpi / 4wpi | 4 days / 4 weeks post injection (COp2) |
| wt / ICE | Control / ICE-treated organoids |
| IN, iN | Inhibitory neurons |
| EN, exN | Excitatory neurons |

## Analysis workflow

### Organoid pipeline (`Organoid_analysis/`)

1. **`01_organoid_QC_processing_integration`** — imports 10x data; MAD-based QC on library size, detected features and mitochondrial content; SCTransform, PCA, clustering, UMAP per dataset; doublet removal (DoubletFinder); Harmony integration; cell-type annotation; removal of stressed-cell clusters.
2. **`02_organoid_Analysis_metacell_Level`** — bootstrapped metacells (pseudobulk) per cell type and dataset group; ICE vs. wt DESeq2 DEA and GO:BP GSEA; timepoint × phenotype interaction testing (LRT + Wald follow-ups) for COp2 4dpi vs. 4wpi.
3. **`03_organoid_Analysis_singleCell_Level`** — single-cell DEA (Seurat `FindMarkers`, Wilcoxon) and GSEA for ICE vs. wt and timepoint comparisons, both across and per cell type.
4. **`04_organoid_ML_predictions_organoids`** — binomial glmnet ICE-vs-wt classifiers per cell type, repeated across random train/test splits, with CV error, ROC/AUC and coefficient-stability diagnostics.

### Human pipeline (`Human_analysis/`)

1. **`01_human_Processing`** — patient/cell filtering and metacell generation for inhibitory (IN) and excitatory (EN) neurons; patient characteristics table.
2. **`02_human_Analysis_metacell_Level`** — variance partitioning by covariates; DESeq2 DEA of Age (adjusting for diagnosis, sample source, sex); GO:BP GSEA.
3. **`02b_human_Analysis_singleCell_Level`** — single-cell negative-binomial GLM (glmGamPoi) for Age, plus GSEA. More sensitive but less conservative than the metacell analysis, since cells are not independent replicates.
4. **`03_human_ML_predictions_human`** — glmnet (gaussian) age regression per cell type, evaluated with (i) repeated 80/20 splits and (ii) leave-patients-out to test generalisation to unseen patients.
5. **`04_human_ML_predictions_human_testOrganoids`** — cross-dataset models: human age model → organoids; organoid ICE classifier → human; human age model → iNeurons.

### Comparison (`compare_organoid_human/`)

**`01_Comparison_human_organoid`** — overlap of glmnet-selected genes between the human and organoid models, concordance of log2 fold changes between human Age DEA and organoid ICE DEA, shared GO terms, and the minimal predictive gene set for organoid models.

**Run order:** organoid 01 → 04 and human 01 → 04 (human `04` needs the organoid metacell outputs), then the comparison notebook.

## Running the analysis

The notebooks are Quarto documents with R code, designed to run on a SLURM cluster.

Edit the `sbatch_*` scripts to select which steps to render (steps are commented in or out). Quarto writes each report next to its notebook; after rendering, move them to `docs/` with:

```bash
bash move_reports_to_docs.sh
```

### Paths and parameters

All input/output paths and analysis thresholds are defined in the `params:` block in each notebook's YAML header (e.g. `helper_functions_path`, `object_storage_path`, `max_age`). They currently point to `/home/vkorob/data/aging_human_organoid_comparisson/...`. If you clone the repository somewhere else, change them in the header.

## Main dependencies

The full conda environment (R, Quarto and all packages with versions) is defined in `environment.yml`:

```bash
conda env create -f environment.yml
```

Main R packages:

- **Single-cell:** Seurat, harmony, DoubletFinder, SingleCellExperiment, scater, glmGamPoi
- **Differential expression:** DESeq2, edgeR, sva
- **Enrichment:** genekitr, geneset, rrvgo
- **Modelling:** glmnet, caret, pROC
- **Data & plotting:** dplyr, tidyr, purrr, stringr, ggplot2, ggpubr, patchwork, pheatmap, EnhancedVolcano, ggVennDiagram, DT, openxlsx

All random steps use `set.seed(1)`; models repeated across splits use a fixed set of seeds.

## Authors

Vladyslav Korobeynyk — Jessberger lab, Brain Research Institute, University of Zurich.
Organoid experiments by Annina Denoth (currently at Roche).

## Citation

Manuscript in preparation.

## License

Copyright (C) 2026 Vladyslav Korobeynyk.

This program is free software: you can redistribute it and/or modify it under the terms of the [GNU General Public License v3.0](LICENSE).
