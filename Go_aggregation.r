# =============================================================================
# GO Pathway Aggregation Analysis — Cross Cell-Line Summary (3D vs 2D)
# Project: AIRC - Cell Line Analysis
# Description: Aggregates GO enrichment results across all cell lines,
#              identifies the most recurrent pathways and visualizes them
#              via a presence/significance heatmap annotated by cancer subtype.
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(openxlsx)
library(ggnewscale)
library(Cairo)

# --- Parameters ---------------------------------------------------------------
cell_features_file <- "data/meta_cell_lines.xlsx"
ontology           <- "BP"    # "BP", "MF", or "CC"
direction          <- "Up"    # "Up" or "Down" — drives both input/output paths and plot titles
base_dir           <- "results"  # parent folder containing GO_analysis_Upgenes and GO_analysis_Downgenes

# Input/output directories set automatically based on direction
input_dir  <- file.path(base_dir, paste0("GO_analysis_", direction, "genes"))
output_dir <- file.path(base_dir, "GO_heatmap")
top_n              <- 20      # number of top recurrent terms to display
term_col           <- "Description"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load cell line metadata --------------------------------------------------
cell_features <- read.xlsx(cell_features_file, sheet = 1)

# --- Load GO results ----------------------------------------------------------
files <- list.files(
  path    = input_dir,
  pattern = paste0("_", ontology, "\\.txt$"),
  full.names = TRUE
)

lista_df <- lapply(files, function(f) {
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
})
names(lista_df) <- gsub(paste0("_", ontology, "\\.txt$"), "", basename(files))

# --- Check empty datasets -----------------------------------------------------
righe_per_df <- sapply(lista_df, nrow)
vuoti <- names(righe_per_df)[righe_per_df == 0]
if (length(vuoti) > 0) {
  message("Empty datasets (will be skipped): ", paste(vuoti, collapse = ", "))
}

# --- Merge all datasets -------------------------------------------------------
df_unico <- bind_rows(lapply(names(lista_df), function(nm) {
  df <- lista_df[[nm]]
  if (nrow(df) == 0) return(NULL)
  df$Dataset <- nm
  df
}))

# --- Count recurrence of each GO term across datasets ------------------------
conteggio_termine <- df_unico %>%
  distinct(Dataset, .data[[term_col]]) %>%
  group_by(.data[[term_col]]) %>%
  summarise(N_dataset = n(), .groups = "drop") %>%
  arrange(desc(N_dataset))

# --- Select top N terms -------------------------------------------------------
top_terms     <- conteggio_termine %>% slice_max(N_dataset, n = top_n)
top_terms     <- top_terms[1:top_n, ]
top_terms_vec <- top_terms[[term_col]]

# --- Build heatmap dataframe --------------------------------------------------
all_datasets <- unique(df_unico$Dataset)

heatmap_df <- df_unico %>%
  filter(.data[[term_col]] %in% top_terms_vec) %>%
  distinct(Dataset, .data[[term_col]], .keep_all = TRUE) %>%
  dplyr::select(Dataset, all_of(term_col), p.adjust) %>%
  complete(
    Dataset = all_datasets,
    !!sym(term_col) := top_terms_vec,
    fill = list(p.adjust = NA)
  ) %>%
  mutate(Significant = case_when(
    is.na(p.adjust)    ~ "Not significant",
    p.adjust < 0.05    ~ "Significant < 0.05 padj",
    TRUE               ~ "Significant < 0.1 padj"
  ))

# Order GO terms by recurrence (most recurrent on top)
term_order <- top_terms %>% arrange(N_dataset) %>% pull(.data[[term_col]])
heatmap_df[[term_col]] <- factor(heatmap_df[[term_col]], levels = term_order)

# Clean dataset names (keep only cell line name before "..")
heatmap_df$Dataset    <- sub("\\.\\..*", "", heatmap_df$Dataset)
heatmap_df$Cell.lines <- heatmap_df$Dataset

# Add cancer subtype from metadata
heatmap_df <- heatmap_df %>%
  left_join(cell_features %>% dplyr::select(Cell.lines, Subtype), by = "Cell.lines")

# Order datasets by subtype
dataset_order <- heatmap_df %>%
  distinct(Dataset, Subtype) %>%
  arrange(Subtype, Dataset) %>%
  pull(Dataset)

heatmap_df$Dataset <- factor(heatmap_df$Dataset, levels = dataset_order)

# Annotation row for subtype strip
annotation_df <- heatmap_df %>%
  distinct(Dataset, Subtype) %>%
  mutate(term_dummy = " ")

# --- Plot heatmap -------------------------------------------------------------
CairoPDF(
  file   = file.path(output_dir, paste0("heatmap_GO_", ontology, "_", direction, "regulated.pdf")),
  width  = 17,
  height = 10,
  family = "Arial"
)

ggplot() +
  
  # Main heatmap — significance fill
  geom_tile(
    data = heatmap_df,
    aes(x = Dataset, y = .data[[term_col]], fill = Significant),
    color = "white", linewidth = 0.5
  ) +
  scale_fill_manual(
    name   = "Significance",
    values = c(
      "Significant < 0.05 padj" = "steelblue",
      "Significant < 0.1 padj"  = "lightblue",
      "Not significant"         = "grey90"
    )
  ) +
  
  # Subtype annotation strip
  new_scale_fill() +
  geom_tile(
    data = annotation_df,
    aes(x = Dataset, y = term_dummy, fill = Subtype),
    height = 0.8
  ) +
  scale_fill_manual(
    name   = "Subtype",
    values = c(
      "LA"  = "#a3a500",
      "H"   = "#f8766d",
      "LB"  = "#00bf7d",
      "TNB" = "#e76bf3",
      "TNA" = "#e79bf7"
    )
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    panel.grid   = element_blank(),
    legend.position = "bottom"
  ) +
  labs(
    title = paste0("Top ", top_n, " GO Terms — ", ontology, " — ", direction, "regulated"),
    x     = "Cell line",
    y     = "GO Term"
  )

dev.off()

message("Heatmap saved: ", file.path(output_dir, paste0("heatmap_GO_", ontology, "_", direction, "regulated.pdf")))
