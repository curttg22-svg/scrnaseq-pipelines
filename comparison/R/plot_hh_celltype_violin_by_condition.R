# =============================================================================
# plot_hh_celltype_violin_by_condition.R
#
# PTCH1, PTCH2, and GLI1 expression by cell type, split by regenerative vs
# non-regenerative condition in the collaborator axolotl dataset.
#
# Because HH genes are sparsely expressed in scRNA-seq (median = 0 in most
# cell types), two complementary views are generated:
#
#   Page 1 per gene — Violin + mean dot, faceted Regen | NonRegen.
#     Shows the distribution shape and mean expression level.
#
#   Final page — Dot plot: % cells expressing (>0) per cell type per condition.
#     The most interpretable view for sparse pathway genes.
#
# Key imbalance: Regenerative Blastema Fibroblast has ~3,800 Regen cells but
# only ~24 NonRegen cells — NonRegen violin will be sparse by construction.
#
# Run from comparison/ directory:
#   source("R/plot_hh_celltype_violin_by_condition.R")
#
# Output: results/hh_celltype_violin_by_condition.pdf
# =============================================================================

library(ggplot2)
library(patchwork)
library(hdf5r)

AXO_COLLAB_DIR <- "data/axo_collab"
RESULTS_DIR    <- "results"

HH_GENES <- c("PTCH1", "PTCH2", "GLI1")

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

COND_COLORS <- c(Regen = "#1565C0", NonRegen = "#B71C1C")

# =============================================================================
# 1. LOAD METADATA + EXPRESSION
# =============================================================================

message("Loading collaborator axolotl metadata...")
meta     <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objmeta.rds"))
gene_idx <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objgene.rds"))[["RNA"]]

meta$condition <- factor(
  ifelse(grepl("^non-regenerating", meta$sample_id), "NonRegen", "Regen"),
  levels = c("Regen", "NonRegen")
)

ct_cond_tab <- table(meta$annotated_celltype, meta$condition)
message("Cell type x condition counts:")
print(ct_cond_tab)

target_rows <- gene_idx[HH_GENES]
target_rows <- target_rows[!is.na(target_rows)]
found_genes <- names(target_rows)
missing     <- setdiff(HH_GENES, found_genes)
if (length(missing) > 0)
  message("Not found in annotation: ", paste(missing, collapse = ", "))
message("Fetching from h5: ", paste(found_genes, collapse = ", "),
        " (rows: ", paste(target_rows, collapse = ", "), ")")

h5   <- H5File$new(file.path(AXO_COLLAB_DIR, "merged_objassay_RNA.h5"), mode = "r")
dset <- h5[["grp/data"]]
axo_expr <- setNames(
  lapply(found_genes, function(g) dset[target_rows[g], ]),
  found_genes
)
h5$close_all()

axo_df <- as.data.frame(axo_expr)
axo_df$cell_type <- meta$annotated_celltype
axo_df$condition <- meta$condition
axo_df <- axo_df[!is.na(axo_df$cell_type), ]

# Order cell types by mean expression of PTCH1 in Regen (descending)
order_gene <- if ("PTCH1" %in% found_genes) "PTCH1" else found_genes[1]
regen_sub  <- axo_df[axo_df$condition == "Regen", ]
ct_order   <- names(sort(tapply(regen_sub[[order_gene]], regen_sub$cell_type, mean),
                         decreasing = TRUE))
axo_df$cell_type <- factor(axo_df$cell_type, levels = ct_order)

# =============================================================================
# 2. VIOLIN PLOTS — one page per gene, faceted by condition
# =============================================================================

make_condition_violin <- function(df, gene) {
  df2  <- df[!is.na(df[[gene]]), ]
  n_re <- sum(df2$condition == "Regen")
  n_nr <- sum(df2$condition == "NonRegen")
  sub  <- paste0("Regen: ", n_re, " cells  |  NonRegen: ", n_nr, " cells  |  ",
                 "mean shown as point; median = 0 for most cell types (sparse gene)")

  ggplot(df2, aes(x = cell_type, y = .data[[gene]], fill = condition)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, color = "grey30") +
    stat_summary(fun = mean, geom = "point", size = 1.8, color = "white",
                 shape = 21, fill = "black", stroke = 0.4) +
    scale_fill_manual(values = COND_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    facet_wrap(~ condition, nrow = 2) +
    labs(
      title    = paste0(gene, " by cell type — Regen vs NonRegen (axolotl collaborator)"),
      subtitle = sub,
      x = NULL, y = paste0(gene, " (log-norm)")
    ) +
    base_theme
}

# =============================================================================
# 3. PERCENTAGE-EXPRESSING DOT PLOT — all genes, all cell types
# =============================================================================

make_pct_dot_plot <- function(df, genes) {
  rows <- do.call(rbind, lapply(genes, function(gene) {
    do.call(rbind, lapply(levels(df$condition), function(cond) {
      sub   <- df[df$condition == cond & !is.na(df[[gene]]), ]
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
  rows$condition <- factor(rows$condition, levels = c("Regen", "NonRegen"))

  ggplot(rows, aes(x = condition, y = cell_type,
                   size = pct_expr, color = mean_expr)) +
    geom_point(alpha = 0.9) +
    scale_size_continuous(name = "% cells\nexpressing", range = c(0.5, 10),
                          breaks = c(1, 5, 10, 20)) +
    scale_color_gradient(low = "#F7F7F7", high = "#B2182B",
                         name = "Mean\nlog-norm") +
    facet_wrap(~ gene, nrow = 1) +
    labs(
      title    = "HH gene expression: % cells expressing and mean level by condition",
      subtitle = paste0("Dot size = % cells with log-norm > 0  |  ",
                        "Color = mean log-norm expression  |  ",
                        "Cell types ordered by mean PTCH1 in Regen"),
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

out_file <- file.path(RESULTS_DIR, "hh_celltype_violin_by_condition.pdf")
pdf(out_file, width = 14, height = 14)

for (gene in found_genes) {
  message("Plotting violin: ", gene)
  print(make_condition_violin(axo_df, gene))
}

message("Plotting summary dot plot...")
print(make_pct_dot_plot(axo_df, found_genes))

dev.off()
message("Saved ", out_file)

# Console summary: % expressing in Regen vs NonRegen per cell type
message("\n--- % cells expressing HH genes by condition ---")
for (gene in found_genes) {
  message("\n", gene, ":")
  for (cond in c("Regen", "NonRegen")) {
    sub   <- axo_df[axo_df$condition == cond, ]
    pct   <- 100 * tapply(sub[[gene]], sub$cell_type, function(x) mean(x > 0))
    pct_s <- sort(pct, decreasing = TRUE)
    msg   <- paste(sprintf("%s=%.1f%%", names(pct_s), pct_s), collapse = ", ")
    message("  ", cond, ": ", msg)
  }
}
