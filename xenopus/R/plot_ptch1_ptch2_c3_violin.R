# =============================================================================
# plot_ptch1_ptch2_c3_violin.R
#
# Violin plots of Ptch1, Ptch2, and C3 expression in Xenopus froglet 7-14dpa
# across keratinocytes, basal epidermis (Epithelial_mixed), and macrophages.
#
# Homeologs are summed per cell before plotting.
# Run from xenopus/ directory: source("R/plot_ptch1_ptch2_c3_violin.R")
#
# Output: results/ptch1_ptch2_c3_violin_froglet_7-14dpa.pdf
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)

results_dir <- "results"

# =============================================================================
# 1. LOAD AND SUBSET
# =============================================================================

message("Loading Xenopus object...")
xen <- readRDS(file.path(results_dir, "xen_merged.rds"))

# Froglet (BL), 7-14dpa only
xen_sub <- subset(xen, condition == "BL" & timepoint == "7-14dpa")

# Cell types of interest
CELL_TYPES <- c("Keratinocyte", "Epithelial_mixed",
                "Macrophage", "Macrophage_2",
                "CT_fibroblast_A", "CT_fibroblast_B", "CT_fibroblast_C",
                "CT_fibroblast_D", "CT_fibroblast_E",
                "CT_progenitor")
xen_sub <- subset(xen_sub, cell_type %in% CELL_TYPES)

# Collapse CT fibroblast subtypes into one group
xen_sub$cell_group <- as.character(xen_sub$cell_type)
xen_sub$cell_group[grepl("^CT_fibroblast", xen_sub$cell_group)] <- "CT_fibroblast"

message("Cells retained: ", ncol(xen_sub))
print(table(xen_sub$cell_type))

# =============================================================================
# 2. SUM HOMEOLOGS PER CELL
# =============================================================================

expr <- GetAssayData(xen_sub, layer = "data")

sum_homeologs <- function(expr_mat, pattern) {
  hits <- grep(pattern, rownames(expr_mat), value = TRUE, perl = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) return(rep(0, ncol(expr_mat)))
  if (length(hits) == 1) return(as.numeric(expr_mat[hits, ]))
  Matrix::colSums(expr_mat[hits, , drop = FALSE])
}

xen_sub$PTCH1 <- sum_homeologs(expr, "(?i)^ptch1[.]")
xen_sub$PTCH2 <- sum_homeologs(expr, "(?i)^ptch2[.]")
xen_sub$GLI1  <- sum_homeologs(expr, "(?i)^gli1[.]")

# =============================================================================
# 3. BUILD PLOT DATA FRAME
# =============================================================================

meta <- xen_sub@meta.data[, c("cell_group", "PTCH1", "PTCH2", "GLI1")]

# Friendly cell type labels
meta$cell_label <- factor(
  meta$cell_group,
  levels = c("Keratinocyte", "Epithelial_mixed",
             "Macrophage", "Macrophage_2",
             "CT_fibroblast", "CT_progenitor"),
  labels = c("Keratinocyte", "Basal epidermis",
             "Macrophage 1", "Macrophage 2",
             "CT fibroblast", "CT progenitor")
)

# Long format
df <- data.frame(
  cell_label = rep(meta$cell_label, 3),
  gene       = rep(c("Ptch1", "Ptch2", "Gli1"), each = nrow(meta)),
  expression = c(meta$PTCH1, meta$PTCH2, meta$GLI1)
)
df <- df[!is.na(df$cell_label), ]
df$gene <- factor(df$gene, levels = c("Ptch1", "Ptch2", "Gli1"))

# =============================================================================
# 4. COLORS
# =============================================================================

CELL_COLORS <- c(
  "Keratinocyte"   = "#6A1B9A",
  "Basal epidermis"= "#AB47BC",
  "Macrophage 1"   = "#E65100",
  "Macrophage 2"   = "#FF8F00",
  "CT fibroblast"  = "#1565C0",
  "CT progenitor"  = "#42A5F5"
)

# =============================================================================
# 5. PLOT
# =============================================================================

p <- ggplot(df, aes(x = cell_label, y = expression, fill = cell_label)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.35, color = "grey30") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               color = "grey20", linewidth = 0.4) +
  scale_fill_manual(values = CELL_COLORS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  labs(
    title    = "Xenopus froglet 7-14 dpa: Ptch1, Ptch2, Gli1 expression",
    subtitle = "Homeologs summed per cell | log-normalized expression",
    x        = NULL,
    y        = "Expression (summed log-norm)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 10),
    strip.text       = element_text(face = "bold.italic", size = 12),
    strip.background = element_rect(fill = "grey93", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 8, color = "grey40")
  )

out_file <- file.path(results_dir, "ptch1_ptch2_c3_violin_froglet_7-14dpa.pdf")
pdf(out_file, width = 10, height = 5)
print(p)
dev.off()
message("Saved ", out_file)
