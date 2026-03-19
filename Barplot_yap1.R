# --- Libraries ----------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(purrr)
library(openxlsx)
library(Cairo)

# --- INPUT GENE AS MARKER PLOT ------------------------------------------------
cells <- c(
  "SUM159PT..2D", "SUM159PT..3D",
  "MX1..2D", "MX1..3D",
  "MDAMB231..2D", "MDAMB231..3D",
  "HCC1395..2D", "HCC1395..3D",
  "CAL51..2D", "CAL51..3D",
  "BT549..2D", "BT549..3D",
  "SUM52PE..2D", "SUM52PE..3D",
  "SUM229PE..2D", "SUM229PE..3D",
  "SUM185PE..2D", "SUM185PE..3D",
  "MFM223..2D", "MFM223..3D",
  "MDAMB468..2D", "MDAMB468..3D",
  "MDAMB436..2D", "MDAMB436..3D",
  "HCC70..2D", "HCC70..3D",
  "HCC1806..2D", "HCC1806..3D",
  "HCC1143..2D", "HCC1143..3D",
  "DU4475..2D", "DU4475..3D",
  "SKBR3..2D", "SKBR3..3D",
  "MDAMB453..2D", "MDAMB453..3D",
  "JIMT1..2D", "JIMT1..3D",
  "AU565..2D", "AU565..3D",
  "MDAMB361..2D", "MDAMB361..3D",
  "BT474..2D", "BT474..3D",
  "ZR751..2D", "ZR751..3D",
  "T47D..2D", "T47D..3D",
  "MDAMB415..2D", "MDAMB415..3D",
  "MCF7..2D", "MCF7..3D",
  "HCC1500..2D", "HCC1500..3D",
  "CAMA1..2D", "CAMA1..3D",
  "BT483..2D", "BT483..3D"
)

# --- IMPORT DATA --------------------------------------------------------------
meta      <- readRDS("path/to/meta_2Dvs3D.rds")
corrected <- readRDS("path/to/corrected_2Dvs3D.rds")
corrected <- corrected@assays@data$corrected / log2(exp(1))
corrected <- corrected[, rownames(meta)]

# --- BINARIZZAZIONE (metodo adattivo per batch) -------------------------------
get_binary <- function(processed, batch) {
  sparsity <- numeric(length(unique(batch)))
  data     <- vector('list', length(unique(batch)))
  
  for (i in seq_along(unique(batch))) {
    data[[i]] <- processed[, batch == unique(batch)[i]]
    sparsity[i] <- sum(data[[i]] > 0) / prod(dim(data[[i]]))
  }
  
  idx <- which.min(sparsity)
  p   <- 1 - sparsity[idx]
  
  for (i in seq_along(data)) {
    if (i != idx) {
      q        <- quantile(data[[i]], probs = p)
      data[[i]] <- data[[i]] > q
    }
  }
  
  data[[idx]] <- data[[idx]] > 0
  
  binary <- do.call(cbind, data)
  binary <- binary[, colnames(processed)]
  
  return(binary)
}

binary <- get_binary(corrected, meta$batch)

# --- CALCOLO PERCENTUALI ------------------------------------------------------
cell_lines <- unique(meta$line_culture)

percent_list <- map(cell_lines, function(cl) {
  col_num    <- which(meta$line_culture == cl)
  df_bin     <- binary[, col_num, drop = FALSE]
  percentage <- (rowSums(df_bin) / ncol(df_bin)) * 100
  return(percentage)
})

# Converti in dataframe e assegna i nomi di colonna
percent_df <- as.data.frame(do.call(cbind, percent_list))
colnames(percent_df) <- cell_lines
percent_df$Gene <- rownames(expr)

# --- SELEZIONE GENE -----------------------------------------------------------
gene_sel <- "YAP1"

# Filtra solo il gene selezionato
mat <- percent_df[rownames(percent_df) == gene_sel, ]

# Filtra e riordina le colonne secondo la lista cells
mat <- mat[, colnames(mat) %in% cells]
mat <- mat[, match(cells, colnames(mat))]

rownames(mat) <- gene_sel
mat$Gene <- NULL

# Trasponi e converti in dataframe
mat     <- t(as.matrix(mat))
df_plot <- data.frame(
  Cell       = rownames(mat),
  Expression = as.numeric(mat[, 1])
)
df_plot$Cell      <- factor(df_plot$Cell, levels = df_plot$Cell)
df_plot$PairColor <- ifelse(seq_len(nrow(df_plot)) %% 2 == 1, "2D", "3D")

# --- AGGIUNGI SUBTYPE ---------------------------------------------------------
cell_features <- xlsx::read.xlsx(
  "/home/tigem/e.aiello/Progetti/AIRC/Data/meta_cell_lines.xlsx",
  sheetIndex = 1
)

# 1. Estrai nome base della cell line
df_plot <- df_plot %>%
  mutate(CellBase = gsub("\\.\\.(2D|3D)$", "", as.character(Cell)))

# 2. Join con cell_features
df_plot <- df_plot %>%
  left_join(cell_features %>% select(Cell.lines, Subtype),
            by = c("CellBase" = "Cell.lines"))

# 3. Ordina per Subtype, poi per CellBase, e aggiorna i livelli del factor
#    *** CHIAVE: levels = unique(as.character(Cell)) allinea i livelli
#    all'ordine delle righe dopo arrange(), evitando sfasamenti nelle strisce ***
df_plot <- df_plot %>%
  arrange(Subtype, CellBase) %>%
  mutate(Cell = factor(Cell, levels = unique(as.character(Cell))))

# 4. Ricalcola pair_positions e pair_labels (DOPO il fix del factor)
pair_positions <- df_plot %>%
  mutate(idx = as.integer(Cell)) %>%
  group_by(CellBase) %>%
  summarise(pos = mean(idx), .groups = "drop") %>%
  arrange(pos) %>%
  pull(pos)

pair_labels <- df_plot %>%
  mutate(idx = as.integer(Cell)) %>%
  group_by(CellBase) %>%
  summarise(pos = mean(idx), .groups = "drop") %>%
  arrange(pos) %>%
  pull(CellBase)

# 5. Striscia Subtype (DOPO il fix del factor)
subtype_blocks <- df_plot %>%
  mutate(idx = as.integer(Cell)) %>%
  group_by(Subtype) %>%
  summarise(xmin = min(idx) - 0.5,
            xmax = max(idx) + 0.5,
            .groups = "drop")

# Palette Subtype
subtype_colors <- c(
  "LA"  = "#a3a500",
  "H"   = "#f8766d",
  "LB"  = "#00bf7d",
  "TNB" = "#e76bf3",
  "TNA" = "#e79bf7"
)

y_stripe      <- -max(df_plot$Expression) * 0.18
stripe_height <-  max(df_plot$Expression) * 0.06

# --- PLOT ---------------------------------------------------------------------
ggplot(df_plot, aes(x = Cell, y = Expression, fill = PairColor)) +
  geom_bar(stat = "identity") +
  
  # Striscia Subtype
  geom_rect(data = subtype_blocks,
            aes(xmin = xmin, xmax = xmax,
                ymin = y_stripe - stripe_height / 2,
                ymax = y_stripe + stripe_height / 2,
                fill = Subtype),
            inherit.aes = FALSE,
            alpha = 0.85) +
  
  # Etichetta Subtype al centro di ogni blocco
  geom_text(data = subtype_blocks %>%
              mutate(xmid = (xmin + xmax) / 2),
            aes(x = xmid, y = y_stripe, label = Subtype),
            inherit.aes = FALSE,
            size = 3, fontface = "bold", color = "white") +
  
  scale_fill_manual(values = c(
    "2D" = "#ADD8E6",
    "3D" = "#1E90FF",
    subtype_colors
  )) +
  
  scale_x_discrete(labels = rep("", nrow(df_plot))) +
  
  annotate("text",
           x     = pair_positions,
           y     = -max(df_plot$Expression) * 0.05,
           label = pair_labels,
           angle = 45, hjust = 1, size = 3.5) +
  
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 5, r = 5, b = 120, l = 5),
    panel.border = element_blank(),
    axis.line    = element_blank()
  ) +
  labs(
    title = paste("Percentage of cells expressing", gene_sel),
    x     = "Cell line",
    y     = "Percentage",
    fill  = "Type"
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))