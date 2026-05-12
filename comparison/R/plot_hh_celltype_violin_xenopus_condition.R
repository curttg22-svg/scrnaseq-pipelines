# =============================================================================
# plot_hh_celltype_violin_xenopus_condition.R
#
# PTCH1, PTCH2, and GLI1 expression by cell type, split by condition in
# Xenopus (Lin 2021):
#
#   BL  = froglet back limb regeneration (regenerative)
#   LBst = tadpole limb bud stages NF50/51/52 (non-regenerative development)
#
# Homeologs summed: PTCH1 = ptch1.L + ptch1.S, etc.
#
# Two views per gene:
#   1. Violin + mean dot, faceted BL | LBst
#   2. Summary dot plot: % cells expressing (>0) by cell type x condition
#
# Run from comparison/ directory:
#   source("R/plot_hh_celltype_violin_xenopus_condition.R")
#
# Output: results/hh_celltype_violin_xenopus_condition.pdf
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)

XEN_RDS     <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
RESULTS_DIR <- "results"

HH_GENES <- c("PTCH1", "PTCH2", "GLI1")

xen_gene_map <- list(
  PTCH1 = c("ptch1.L", "ptch1.S"),
  PTCH2 = c("ptch2.L", "ptch2.S"),
  GLI1  = c("gli1.L",  "gli1.S")
)

# =============================================================================
# THEME
# =============================================================================

base_theme <- theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey93", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, color = "grey40"),
    legend.position    = "none"
  )

COND_COLORS <- c(BL = "#1565C0", LBst = "#2E7D32")

# =============================================================================
# HELPER — extract gene matrix from Seurat v5 split layers
# =============================================================================

fetch_genes_split <- function(seu, genes, meta_cols = character(0)) {
  all_layers  <- Layers(seu)
  data_layers <- all_layers[startsWith(all_layers, "data.") | all_layers == "data"]
  expr_list <- lapply(data_layers, function(ln) {
    mat  <- LayerData(seu, layer = ln)
    hits <- intersect(genes, rownames(mat))
    if (length(hits) == 0) return(NULL)
    df <- as.data.frame(t(as.matrix(mat[hits, , drop = FALSE])))
    for (g in setdiff(genes, hits)) df[[g]] <- 0
    df[, genes, drop = FALSE]
  })
  expr_df <- do.call(rbind, Filter(Negate(is.null), expr_list))
  if (length(meta_cols) > 0) {
    meta    <- seu@meta.data[rownames(expr_df), meta_cols, drop = FALSE]
    expr_df <- cbind(expr_df, meta)
  }
  expr_df
}

# =============================================================================
# 1. LOAD XENOPUS + EXTRACT HH GENES
# =============================================================================

message("Loading Xenopus RDS...")
xen <- readRDS(path.expand(XEN_RDS))

fetch_genes <- unique(unlist(xen_gene_map))
message("Fetching homeolog pairs: ", paste(fetch_genes, collapse = ", "))

xen_raw <- fetch_genes_split(xen, fetch_genes, meta_cols = c("condition", "cell_type"))
rm(xen); gc()

# Sum homeologs
xen_df <- data.frame(
  condition = xen_raw$condition,
  cell_type = xen_raw$cell_type,
  stringsAsFactors = FALSE
)
for (gene in HH_GENES) {
  homeologs        <- intersect(xen_gene_map[[gene]], colnames(xen_raw))
  xen_df[[gene]]   <- if (length(homeologs) > 0)
    rowSums(xen_raw[, homeologs, drop = FALSE]) else 0
}

xen_df$condition <- factor(xen_df$condition, levels = c("BL", "LBst"))
xen_df <- xen_df[!is.na(xen_df$cell_type) & !is.na(xen_df$condition), ]

ct_cond_tab <- table(xen_df$cell_type, xen_df$condition)
message("\nCell type x condition counts:")
print(ct_cond_tab)

# Order cell types by mean PTCH1 in BL (descending)
bl_sub   <- xen_df[xen_df$condition == "BL", ]
ct_order <- names(sort(tapply(bl_sub[["PTCH1"]], bl_sub$cell_type, mean),
                       decreasing = TRUE))
xen_df$cell_type <- factor(xen_df$cell_type, levels = ct_order)

# =============================================================================
# 2. VIOLIN PLOTS — one page per gene
# =============================================================================

make_condition_violin <- function(df, gene) {
  df2  <- df[!is.na(df[[gene]]), ]
  n_bl   <- sum(df2$condition == "BL",   na.rm = TRUE)
  n_lb   <- sum(df2$condition == "LBst", na.rm = TRUE)
  sub    <- paste0("BL (regenerative froglet): ", n_bl, " cells  |  ",
                   "LBst (non-regen tadpole): ", n_lb, " cells  |  ",
                   "mean shown as point")

  ggplot(df2, aes(x = cell_type, y = .data[[gene]], fill = condition)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, color = "grey30") +
    stat_summary(fun = mean, geom = "point", size = 1.8,
                 shape = 21, fill = "black", color = "white", stroke = 0.4) +
    scale_fill_manual(values = COND_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    facet_wrap(~ condition, nrow = 2) +
    labs(
      title    = paste0(gene, " by cell type — BL vs LBst (Xenopus, Lin 2021)"),
      subtitle = sub,
      x = NULL, y = paste0(gene, " (log-norm, .L+.S summed)")
    ) +
    base_theme
}

# =============================================================================
# 3. SUMMARY DOT PLOT — % expressing
# =============================================================================

make_pct_dot_plot <- function(df, genes, ct_order) {
  rows <- do.call(rbind, lapply(genes, function(gene) {
    do.call(rbind, lapply(levels(df$condition), function(cond) {
      sub   <- df[df$condition == cond & !is.na(df[[gene]]), ]
      if (nrow(sub) == 0) return(NULL)
      n_tot <- tapply(sub[[gene]], sub$cell_type, length)
      n_pos <- tapply(sub[[gene]], sub$cell_type, function(x) sum(x > 0))
      pct   <- 100 * n_pos / n_tot
      mean_expr <- tapply(sub[[gene]], sub$cell_type, mean)
      data.frame(
        gene      = gene,
        condition = cond,
        cell_type = names(pct),
        pct_expr  = as.numeric(pct),
        mean_expr = as.numeric(mean_expr),
        n_cells   = as.integer(n_tot),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rows$cell_type <- factor(rows$cell_type, levels = rev(ct_order))
  rows$gene      <- factor(rows$gene, levels = genes)
  rows$condition <- factor(rows$condition, levels = c("BL", "LBst"))

  ggplot(rows, aes(x = condition, y = cell_type,
                   size = pct_expr, color = mean_expr)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(name = "% cells\nexpressing", range = c(0.5, 10),
                          breaks = c(1, 5, 10, 20)) +
    scale_color_gradient(low = "#F7F7F7", high = "#B2182B",
                         name = "Mean\nlog-norm") +
    facet_wrap(~ gene, nrow = 1) +
    labs(
      title    = "HH gene expression: % cells expressing and mean level by condition (Xenopus)",
      subtitle = "BL = regenerative froglet back limb  |  LBst = non-regenerative tadpole limb bud (NF50-52)\nDot size = % cells with log-norm > 0  |  Color = mean log-norm  |  Cell types ordered by mean PTCH1 in BL",
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x      = element_text(size = 9),
      axis.text.y      = element_text(size = 8),
      strip.text       = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "grey93", color = NA),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 8, color = "grey40")
    )
}

# =============================================================================
# 4. SAVE PDF
# =============================================================================

out_file <- file.path(RESULTS_DIR, "hh_celltype_violin_xenopus_condition.pdf")
pdf(out_file, width = 14, height = 14)

for (gene in HH_GENES) {
  message("Plotting violin: ", gene)
  print(make_condition_violin(xen_df, gene))
}

message("Plotting summary dot plot...")
print(make_pct_dot_plot(xen_df, HH_GENES, ct_order))

dev.off()
message("Saved ", out_file)

# Console summary
message("\n--- % cells expressing HH genes by condition (Xenopus) ---")
for (gene in HH_GENES) {
  message("\n", gene, " (homeologs summed):")
  for (cond in c("BL", "LBst")) {
    sub  <- xen_df[xen_df$condition == cond, ]
    pct  <- 100 * tapply(sub[[gene]], sub$cell_type, function(x) mean(x > 0))
    pcts <- sort(pct, decreasing = TRUE)
    msg  <- paste(sprintf("%s=%.1f%%", names(pcts), pcts), collapse = ", ")
    message("  ", cond, ": ", msg)
  }
}
