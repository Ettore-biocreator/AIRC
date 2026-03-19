# AIRC — Functional Analysis of 3D vs 2D Breast Cancer Cell Lines

> GO enrichment analysis pipeline for scRNA-seq data from 29 breast cancer cell lines grown in 2D and 3D (organoids). Identifies recurrent pathways across cell lines and visualizes single-gene expression profiles annotated by cancer subtype.

---

## Project context

This repository contains the functional analysis component of a larger project studying transcriptional differences between 2D and 3D (organoid) culture systems across 29 breast cancer cell lines using single-cell RNA-seq data.

Differentially expressed genes (DEGs) were obtained from a linear mixed model (`lme4`). Because the model flattens logFC values, the top 1% most up- or downregulated genes per cell line were selected as input for the functional analysis.

---

## Repository structure

```
AIRC/
├── GO_analysis.R                      # GO enrichment analysis (up and downregulated)
├── GO_pathway_aggregation.R           # Cross-dataset pathway recurrence and heatmap
├── YAP1_expression_barplot.R          # Single-gene expression visualization
└── README.md
```

> **Note:** Input data files are not included in this repository. Each script contains a parameters section at the top where input file paths can be set.

---

## Scripts

### 1. `GO_analysis.R`
Performs GO enrichment analysis on the top 1% up- or downregulated DEGs for each cell line.

**Key parameters:**
| Parameter | Description |
|---|---|
| `direction` | `"Up"` or `"Down"` — selects gene direction and output folder |
| `input_file` | Path to the Excel file with DEG results |
| `padj_cutoff` | Adjusted p-value threshold for DEG filtering (default: `0.05`) |
| `pct_genes` | Fraction of genes to select per cell line (default: `0.01`) |

**Output:** One `.txt` file per cell line per ontology (`BP`, `MF`, `CC`) in `results/GO_analysis_Upgenes/` or `results/GO_analysis_Downgenes/`.

---

### 2. `GO_pathway_aggregation.R`
Aggregates GO results across all cell lines, identifies the most recurrent pathways and produces a significance heatmap annotated by cancer subtype.

**Key parameters:**
| Parameter | Description |
|---|---|
| `direction` | `"Up"` or `"Down"` — drives input folder, output filename and plot title |
| `base_dir` | Parent folder containing `GO_analysis_Upgenes/` and `GO_analysis_Downgenes/` |
| `ontology` | GO ontology to analyze: `"BP"`, `"MF"`, or `"CC"` |
| `top_n` | Number of top recurrent terms to display (default: `20`) |

**Output:** PDF heatmap in `results/GO_heatmap/`.

---

### 3. `YAP1_expression_barplot.R`
Computes the percentage of cells expressing a gene of interest per cell line (2D and 3D), using a batch-aware binarization strategy, and visualizes the result as a barplot annotated by cancer subtype.

**Key parameters:**
| Parameter | Description |
|---|---|
| `gene_sel` | Gene of interest (default: `"YAP1"`) |
| `input_rds` | Path to corrected expression matrix (`.rds`) |
| `meta_rds` | Path to single-cell metadata (`.rds`) |

**Binarization:** Expression is binarized using a batch-aware threshold — the sparsest batch sets the reference percentile, which is then applied to all other batches. This approach is more robust than a fixed threshold.

**Output:** PDF barplot in `results/YAP1_expression/`.

---

## Dependencies

```r
install.packages(c("dplyr", "ggplot2", "tidyr", "openxlsx", "purrr", "Cairo"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "rrvgo"))

install.packages("ggnewscale")
```

---

## Notes

- All scripts use **relative paths** — set your working directory to the root of the repository before running.
- Output directories are created automatically if they do not exist.
- To avoid namespace conflicts between packages, `dplyr::select()` and `dplyr::filter()` are called explicitly throughout. It is recommended to add `library(conflicted)` at the top of each session.
