# =============================================================================
# plot_housekeeping_violin.R
#
# Violin plots of ACTB and GAPDH expression across all three species datasets
# to validate normalization consistency across timepoints.
#
# Reads from RDS objects directly; extracts only the 2-3 target genes from
# split layers without calling JoinLayers on the full object.
#
# Run from comparison/ directory: source("R/plot_housekeeping_violin.R")
#
# Output: results/housekeeping_violin.pdf
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)

AXO_RDS   <- "~/Desktop/Axolotl_scrnaseq_project_V3/results/axo_integrated.rds"
XEN_RDS   <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
MOUSE_RDS <- "~/Desktop/scrnaseq-pipelines/mouse-digit/results/mouse_merged.rds"
RESULTS_DIR <- "results"

# =============================================================================
# HELPERS
# =============================================================================

# Efficiently pull a few genes from a Seurat v5 object with split data layers.
# Returns a data.frame with one row per cell, columns = gene names (plus
# any metadata columns requested).
fetch_genes_split <- function(seu, genes, meta_cols = character(0)) {
  all_layers <- Layers(seu)
  data_layers <- all_layers[startsWith(all_layers, "data.") | all_layers == "data"]

  expr_list <- lapply(data_layers, function(ln) {
    mat   <- LayerData(seu, layer = ln)
    hits  <- intersect(genes, rownames(mat))
    if (length(hits) == 0) return(NULL)
    df <- as.data.frame(t(as.matrix(mat[hits, , drop = FALSE])))
    missing <- setdiff(genes, hits)
    for (g in missing) df[[g]] <- 0
    df[, genes, drop = FALSE]
  })

  expr_df <- do.call(rbind, Filter(Negate(is.null), expr_list))

  if (length(meta_cols) > 0) {
    meta <- seu@meta.data[rownames(expr_df), meta_cols, drop = FALSE]
    expr_df <- cbind(expr_df, meta)
  }
  expr_df
}

sum_homeologs_from_mat <- function(mat, pattern) {
  hits <- grep(pattern, rownames(mat), value = TRUE, perl = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) return(rep(0, ncol(mat)))
  if (length(hits) == 1) return(as.numeric(mat[hits, ]))
  Matrix::colSums(mat[hits, , drop = FALSE])
}

violin_theme <- theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 9),
    strip.text       = element_text(face = "bold.italic", size = 11),
    strip.background = element_rect(fill = "grey93", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, color = "grey40")
  )

make_tp_violin <- function(df, value_col, tp_col, gene_label, title, subtitle,
                            tp_levels, tp_colors) {
  df2 <- df[!is.na(df[[tp_col]]), ]
  df2[[tp_col]] <- factor(df2[[tp_col]], levels = tp_levels)
  ggplot(df2, aes(x = .data[[tp_col]], y = .data[[value_col]],
                  fill = .data[[tp_col]])) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.35, color = "grey30") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 color = "grey20", linewidth = 0.4) +
    scale_fill_manual(values = tp_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = title, subtitle = subtitle,
         x = NULL, y = paste0(gene_label, " (log-norm)")) +
    violin_theme
}

# =============================================================================
# 1. AXOLOTL
# =============================================================================

message("Loading Axolotl...")
axo <- readRDS(path.expand(AXO_RDS))
axo$paper_cluster <- trimws(axo$paper_cluster)

AXO_TP  <- c("Intact", "WH", "EB", "MB")
AXO_COL <- c(Intact = "#1A237E", WH = "#283593", EB = "#1565C0", MB = "#42A5F5")

axo_df_raw <- fetch_genes_split(axo, "ACTB", meta_cols = "timepoint")
rm(axo); gc()

axo_df <- data.frame(
  timepoint = factor(axo_df_raw$timepoint, levels = AXO_TP),
  ACTB      = axo_df_raw[["ACTB"]]
)

n_axo <- table(droplevels(axo_df$timepoint))
p_axo_actb <- make_tp_violin(axo_df, "ACTB", "timepoint",
  gene_label = "ACTB",
  title    = "Axolotl (Leigh et al. 2018)",
  subtitle = paste0("n cells per timepoint: ",
                    paste(paste0(names(n_axo), "=", n_axo), collapse = ", "),
                    "  |  GAPDH: not annotated in Trinity assembly"),
  tp_levels = AXO_TP, tp_colors = AXO_COL)

# =============================================================================
# 2. XENOPUS
# =============================================================================

message("Loading Xenopus...")
xen <- readRDS(path.expand(XEN_RDS))

XEN_BL_TP   <- c("0dpa","3dpa","7-14dpa","14dpa","14-52dpa")
XEN_LBST_TP <- c("NF50","NF51","NF52")
XEN_BL_COL  <- c("0dpa"="#E65100","3dpa"="#F57C00","7-14dpa"="#FB8C00",
                  "14dpa"="#FFA726","14-52dpa"="#FFCC80")
XEN_LB_COL  <- c(NF50="#1B5E20", NF51="#388E3C", NF52="#66BB6A")

# Pull just the homeolog genes efficiently from split layers
xen_genes <- c("actb.L", "gapdh.L", "gapdh.S")
xen_df_raw <- fetch_genes_split(xen, xen_genes,
                                 meta_cols = c("condition", "timepoint"))
rm(xen); gc()

xen_df_raw$ACTB  <- xen_df_raw[["actb.L"]]
xen_df_raw$GAPDH <- xen_df_raw[["gapdh.L"]] + xen_df_raw[["gapdh.S"]]

# BL (froglet) subset
bl_df <- xen_df_raw[xen_df_raw$condition == "BL", ]
n_bl  <- table(factor(bl_df$timepoint, levels = XEN_BL_TP))

p_bl_actb <- make_tp_violin(bl_df, "ACTB", "timepoint", "ACTB",
  title    = "Xenopus froglet - BL (Lin et al. 2021)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_bl), "=", n_bl), collapse=", "),
                    "  |  actb.L only (actb.S not detected)"),
  tp_levels = XEN_BL_TP, tp_colors = XEN_BL_COL)

p_bl_gapdh <- make_tp_violin(bl_df, "GAPDH", "timepoint", "GAPDH (summed .L+.S)",
  title    = "Xenopus froglet - BL (Lin et al. 2021)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_bl), "=", n_bl), collapse=", ")),
  tp_levels = XEN_BL_TP, tp_colors = XEN_BL_COL)

# LBst (tadpole) subset
lb_df <- xen_df_raw[xen_df_raw$condition == "LBst", ]
n_lb  <- table(factor(lb_df$timepoint, levels = XEN_LBST_TP))

p_lb_actb <- make_tp_violin(lb_df, "ACTB", "timepoint", "ACTB",
  title    = "Xenopus tadpole - LBst (Lin et al. 2021)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_lb), "=", n_lb), collapse=", ")),
  tp_levels = XEN_LBST_TP, tp_colors = XEN_LB_COL)

p_lb_gapdh <- make_tp_violin(lb_df, "GAPDH", "timepoint", "GAPDH (summed .L+.S)",
  title    = "Xenopus tadpole - LBst (Lin et al. 2021)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_lb), "=", n_lb), collapse=", ")),
  tp_levels = XEN_LBST_TP, tp_colors = XEN_LB_COL)

rm(xen_df_raw, bl_df, lb_df); gc()

# =============================================================================
# 3. MOUSE
# =============================================================================

message("Loading Mouse...")
mus <- readRDS(path.expand(MOUSE_RDS))

MUS_TP  <- c("0dpa","11dpa","12dpa","14dpa","17dpa")
MUS_COL <- c("0dpa"="#1B5E20","11dpa"="#2E7D32","12dpa"="#388E3C",
              "14dpa"="#66BB6A","17dpa"="#A5D6A7")

mus_df_raw <- fetch_genes_split(mus, c("Actb","Gapdh"),
                                 meta_cols = "timepoint")
rm(mus); gc()

n_mus <- table(factor(mus_df_raw$timepoint, levels = MUS_TP))

p_mus_actb <- make_tp_violin(mus_df_raw, "Actb", "timepoint", "ACTB",
  title    = "Mouse digit tip (Johnson, Masias & Lehoczky 2020)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_mus),"=",n_mus), collapse=", ")),
  tp_levels = MUS_TP, tp_colors = MUS_COL)

p_mus_gapdh <- make_tp_violin(mus_df_raw, "Gapdh", "timepoint", "GAPDH",
  title    = "Mouse digit tip (Johnson, Masias & Lehoczky 2020)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_mus),"=",n_mus), collapse=", ")),
  tp_levels = MUS_TP, tp_colors = MUS_COL)

rm(mus_df_raw); gc()

# =============================================================================
# 4. COMBINE AND SAVE
# =============================================================================

# Page 1: ACTB across all datasets
actb_page <- (p_axo_actb / p_bl_actb / p_lb_actb / p_mus_actb) +
  plot_annotation(
    title    = "ACTB expression by timepoint - normalization QC",
    subtitle = "Stable housekeeping gene expression validates log-normalization across timepoints",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

# Page 2: GAPDH (Xenopus + Mouse only - absent from axolotl Trinity assembly)
gapdh_page <- (p_bl_gapdh / p_lb_gapdh / p_mus_gapdh) +
  plot_annotation(
    title    = "GAPDH expression by timepoint - normalization QC",
    subtitle = "GAPDH not annotated in axolotl Trinity assembly (Leigh et al. 2018)",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

out_file <- file.path(RESULTS_DIR, "housekeeping_violin.pdf")
pdf(out_file, width = 10, height = 16)
print(actb_page)
print(gapdh_page)
dev.off()
message("Saved ", out_file)
