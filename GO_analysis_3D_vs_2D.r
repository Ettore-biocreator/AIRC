# =============================================================================
# GO Enrichment Analysis - Upregulated and Downregulated Genes (3D vs 2D)
# Project: AIRC - Cell Line Analysis
# Description: For each cell line, selects the top 1% up- or downregulated DEGs
#              (from a linear mixed model) and performs GO enrichment analysis
#              across BP, MF, and CC ontologies.
# Note: 1% threshold is used because logFC values are flattened by the lme4 model
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)

# --- Parameters ---------------------------------------------------------------
input_file  <- "data/lmtren3Dvs2D_lineByLine.xlsx"
padj_cutoff <- 0.05
pct_genes   <- 0.01  # top 1% most up/downregulated genes per cell line

# Set direction: "up" for upregulated, "down" for downregulated
direction  <- "down"

# Output directory is set automatically based on direction
output_dir <- ifelse(
  direction == "up",
  "results/GO_analysis_Upgenes",
  "results/GO_analysis_Downgenes"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data and select genes -----------------------------------------------
sheets  <- getSheetNames(input_file)
n_sheet <- length(sheets)

list_gene <- list()

for (i in seq_len(n_sheet)) {
  dataset    <- read.xlsx(input_file, sheet = i)
  sheet_name <- sheets[i]
  
  # Filter by adjusted p-value
  dataset <- dataset[dataset$adjpvalue < padj_cutoff, ]
  
  # Order by logFC: ascending for down (most negative first), descending for up
  decreasing <- ifelse(direction == "up", TRUE, FALSE)
  dataset    <- dataset[order(dataset$logFC, decreasing = decreasing), ]
  
  # Select top 1% genes
  n_gene <- max(1, round(nrow(dataset) * pct_genes))
  list_gene[[sheet_name]] <- dataset$names[1:n_gene]
}

# --- Plot: number of selected genes per cell line -----------------------------
gene_sizes <- sapply(list_gene, length)

df_sizes <- data.frame(
  Dataset = names(gene_sizes),
  N_Genes = as.numeric(gene_sizes)
)

plot_title <- ifelse(
  direction == "up",
  "Number of selected upregulated genes per cell line",
  "Number of selected downregulated genes per cell line"
)

ggplot(df_sizes, aes(x = Dataset, y = N_Genes)) +
  geom_bar(stat = "identity", fill = ifelse(direction == "up", "tomato", "steelblue")) +
  geom_text(aes(label = N_Genes), vjust = -0.5, size = 3.5) +
  theme_minimal() +
  labs(title = plot_title, x = "Cell line", y = "Number of genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "n_genes_per_cellline.pdf"), width = 10, height = 5)

# --- GO enrichment analysis ---------------------------------------------------
for (i in seq_along(list_gene)) {
  name_cell <- names(list_gene)[i]
  genes     <- unlist(list_gene[i])
  
  for (ontology in c("BP", "MF", "CC")) {
    enrich_result <- enrichGO(
      gene          = genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = ontology,
      pvalueCutoff  = 0.1,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.2,
      minGSSize     = 10,
      maxGSSize     = 500,
      readable      = FALSE,
      pool          = FALSE
    )
    
    out_file <- file.path(output_dir, paste0(name_cell, "_", ontology, ".txt"))
    write.table(as.data.frame(enrich_result), file = out_file, sep = "\t", quote = FALSE)
  }
  
  message("Done: ", name_cell)
}
