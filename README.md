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
├── environment.yml                       # conda: R 4.3.3, Quarto, renv, build tools
├── renv.lock                             # exact versions of all R packages (renv)
├── renv/                                 # renv activation script & settings (package library not tracked)
├── .Rprofile                             # activates renv (also one in each analysis folder)
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
│   └── organoid_analysis_output/         # objects, DEA, GSEA and glmnet results (not tracked)
├── Human_analysis/
│   ├── 01_human_Processing.qmd
│   ├── 02_human_Analysis_metacell_Level.qmd
│   ├── 03_human_ML_predictions_human.qmd
│   ├── 04_human_ML_predictions_human_testOrganoids.qmd
│   ├── 01_human_extended_Processing.qmd
│   ├── 02_human_extended_Analysis_pseudobulk_Level.qmd
│   ├── 03_human_extended_ML_predictions_human.qmd
│   ├── 04_human_extended_ML_predictions_human_testOrganoids.qmd
│   ├── 05_human_cross_cohort_predictions.qmd
│   ├── sbatch_analysis_human
│   ├── logs/                             # SLURM logs (contents not tracked)
│   ├── human_analysis_output/            # original cohort: objects, DEA, GSEA, glmnet results (not tracked)
│   └── human_extended_analysis_output/   # extended cohort: objects, DEA, GSEA, glmnet results (not tracked)
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
| Human (extended) | snRNA-seq of human DLPFC across the lifespan (~280 donors), provided as donor-level pseudobulk counts per cell subtype; all excitatory and all inhibitory neuron subtypes are summed per donor, donors with ≤ 100 cells of a class are removed, and analyses use donors aged 30–80 | [Synapse syn52396927](https://www.synapse.org/Synapse:syn52396927) (`lifespan_pseudobulk.rds`) |
| iNeurons | Bulk RNA-seq of induced neurons directly reprogrammed from fibroblasts of donors of different ages (22 runs, ERR668328–ERR668411), used as an additional test set for the human age model. FASTQs were processed to a count table with [ARMOR](https://github.com/csoneson/ARMOR) ([Orjuela et al., 2019](https://doi.org/10.1534/g3.119.400185)) and stored as `edgeR_dge_induced_Neurons.rds` | Paper: [Mertens et al., Cell Stem Cell (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5929130/) |

Raw data, `.rds` objects and job logs are **not tracked** in git (see `.gitignore`). They are expected inside the project folder (`organoid_data/`, `human_data/`, `human_data_extended/`, etc.).

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

**Original cohort** (Nature Genetics 2024)

1. **`01_human_Processing`** — patient/cell filtering and metacell generation for inhibitory (IN) and excitatory (EN) neurons; patient characteristics table.
2. **`02_human_Analysis_metacell_Level`** — variance partitioning by covariates; DESeq2 DEA of Age (adjusting for diagnosis, sample source, sex); GO:BP GSEA.
3. **`03_human_ML_predictions_human`** — glmnet (gaussian) age regression per cell type, evaluated with (i) repeated 80/20 splits and (ii) leave-patients-out to test generalisation to unseen patients.
4. **`04_human_ML_predictions_human_testOrganoids`** — cross-dataset models: human age model → organoids; organoid ICE classifier → human; human age model → iNeurons.

**Extended lifespan cohort** (Synapse syn52396927)

1. **`01_human_extended_Processing`** — sums all EN and all IN subtype pseudobulks per donor, records the number of cells per donor, removes donors with few cells and LINC genes, and saves one `SingleCellExperiment` per neuron class.
2. **`02_human_extended_Analysis_pseudobulk_Level`** — variance partitioning; donor-level DESeq2 DEA of Age (adjusting for sex, ethnicity, brain bank, control status and log number of cells); GO:BP GSEA.
3. **`03_human_extended_ML_predictions_human`** — glmnet (gaussian) age regression per cell type, evaluated with repeated 80/20 splits over donors (one sample per donor, so every test donor is unseen).
4. **`04_human_extended_ML_predictions_human_testOrganoids`** — the cross-dataset models of the original `04`, with the extended cohort as the human dataset.
5. **`05_human_cross_cohort_predictions`** — trains an age model on one human cohort and tests it on the other (both directions), compared with the within-cohort cross-validation error; both cohorts as donor-level pseudobulks.

### Comparison (`compare_organoid_human/`)

**`01_Comparison_human_organoid`** — overlap of glmnet-selected genes between the human and organoid models, concordance of log2 fold changes between human Age DEA and organoid ICE DEA, shared GO terms, and the minimal predictive gene set for organoid models.

**Run order:** organoid 01 → 04, human original 01 → 04 and human extended 01 → 05 (both human `04` notebooks need the organoid metacell outputs; `05` needs both human cohorts), then the comparison notebook.

## Running the analysis

The notebooks are Quarto documents with R code, designed to run on a SLURM cluster. Activate the environment first (see [Software environment](#software-environment)); renv is then picked up automatically when rendering from within the repository.

Edit the `sbatch_*` scripts to select which steps to render (steps are commented in or out). Quarto writes each report next to its notebook; after rendering, move them to `docs/` with:

```bash
bash move_reports_to_docs.sh
```

### Paths and parameters

All input/output paths and analysis thresholds are defined in the `params:` block in each notebook's YAML header (e.g. `helper_functions_path`, `object_storage_path`, `max_age`). They currently point to `/home/vkorob/data/aging_human_organoid_comparisson/...`. If you clone the repository somewhere else, change them in the header.

## Software environment

R packages are pinned with [renv](https://rstudio.github.io/renv/): `renv.lock` records the exact version and source (CRAN, Bioconductor 3.18, GitHub) of every package used. R itself, Quarto and the build tools come from a small conda environment (`environment.yml`).

To recreate the environment:

```bash
git clone https://github.com/vkorobeynyk/aging_human_organoid_comparisson.git
cd aging_human_organoid_comparisson

conda env create -f environment.yml          # R 4.3.3, Quarto, renv, compilers
conda activate aging_human_organoid_comparisson

Rscript -e 'renv::restore()'                 # installs all R packages at the locked versions
```

The first `renv::restore()` compiles many packages and can take a while. Afterwards every R session and every Quarto render inside the repository uses the project library automatically: the `.Rprofile` in the repository root and in each analysis folder activates renv.

Main R packages:

- **Single-cell:** Seurat, SCTransform (v2 regularisation, `nb_offset` method), harmony, DoubletFinder, SingleCellExperiment, scater
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
