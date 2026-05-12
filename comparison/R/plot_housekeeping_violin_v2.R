# =============================================================================
# plot_housekeeping_violin_v2.R
#
# QC violin plots of housekeeping gene expression across conditions/timepoints.
# Validates normalization consistency across the compendium (Leigh removed;
# collaborator axolotl replaces it).
#
# Collaborator axolotl: ACTG1 (ACTB absent from annotation) + GAPDH
#   across Regen (dpa3/dpa14/dpa23) and NonRegen (limb) groups.
# Xenopus (Lin 2021): actb.L + gapdh.L+.S across BL dpa timepoints and
#   LBst NF stages.
# Mouse digit tip: Actb + Gapdh across dpa timepoints.
#
# Run from comparison/ directory: source("R/plot_housekeeping_violin_v2.R")
#
# Output: results/housekeeping_violin_v2.pdf
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(hdf5r)

AXO_COLLAB_DIR <- "data/axo_collab"
XEN_RDS        <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
MOUSE_RDS      <- "~/Desktop/scrnaseq-pipelines/mouse-digit/results/mouse_merged.rds"
RESULTS_DIR    <- "results"

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
# 1. COLLABORATOR AXOLOTL — read from ShinyCell h5
# =============================================================================

message("Loading collaborator axolotl...")
meta     <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objmeta.rds"))
gene_idx <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objgene.rds"))[["RNA"]]

# ACTB absent; use ACTG1 (gamma-actin) + GAPDH
target_genes <- c(ACTG1 = "ACTG1", GAPDH = "GAPDH")
target_rows  <- gene_idx[target_genes]
message("  ACTG1 row: ", target_rows["ACTG1"],
        "  GAPDH row: ", target_rows["GAPDH"])

h5   <- H5File$new(file.path(AXO_COLLAB_DIR, "merged_objassay_RNA.h5"), mode = "r")
dset <- h5[["grp/data"]]
axo_expr <- list(
  ACTG1 = dset[target_rows["ACTG1"], ],
  GAPDH = dset[target_rows["GAPDH"], ]
)
h5$close_all()

meta$condition <- ifelse(grepl("^non-regenerating", meta$sample_id),
                         "NonRegen", "Regen")
meta$group <- paste0(meta$condition, "_", meta$timepoint)

axo_df <- data.frame(
  group = factor(meta$group,
                 levels = c("Regen_dpa3", "Regen_dpa14", "Regen_dpa23",
                            "NonRegen_limb")),
  ACTG1 = axo_expr$ACTG1,
  GAPDH = axo_expr$GAPDH
)

AXO_TP  <- c("Regen_dpa3", "Regen_dpa14", "Regen_dpa23", "NonRegen_limb")
AXO_COL <- c(Regen_dpa3 = "#1565C0", Regen_dpa14 = "#1E88E5",
              Regen_dpa23 = "#42A5F5", NonRegen_limb = "#B71C1C")

n_axo <- table(droplevels(axo_df$group))

p_axo_actg1 <- make_tp_violin(axo_df, "ACTG1", "group", "ACTG1",
  title    = "Axolotl (collaborator)",
  subtitle = paste0("n cells per group: ",
                    paste(paste0(names(n_axo), "=", n_axo), collapse = ", "),
                    "  |  ACTB not annotated; ACTG1 used as actin proxy"),
  tp_levels = AXO_TP, tp_colors = AXO_COL)

p_axo_gapdh <- make_tp_violin(axo_df, "GAPDH", "group", "GAPDH",
  title    = "Axolotl (collaborator)",
  subtitle = paste0("n cells per group: ",
                    paste(paste0(names(n_axo), "=", n_axo), collapse = ", ")),
  tp_levels = AXO_TP, tp_colors = AXO_COL)

rm(axo_expr, axo_df); gc()

# =============================================================================
# 2. XENOPUS
# =============================================================================

message("Loading Xenopus...")
xen <- readRDS(path.expand(XEN_RDS))

XEN_BL_TP  <- c("0dpa","3dpa","7-14dpa","14dpa","14-52dpa")
XEN_LB_TP  <- c("NF50","NF51","NF52")
XEN_BL_COL <- c("0dpa"="#E65100","3dpa"="#F57C00","7-14dpa"="#FB8C00",
                 "14dpa"="#FFA726","14-52dpa"="#FFCC80")
XEN_LB_COL <- c(NF50="#1B5E20", NF51="#388E3C", NF52="#66BB6A")

xen_raw <- fetch_genes_split(xen, c("actb.L","gapdh.L","gapdh.S"),
                              meta_cols = c("condition","timepoint"))
rm(xen); gc()

xen_raw$ACTB  <- xen_raw[["actb.L"]]
xen_raw$GAPDH <- xen_raw[["gapdh.L"]] + xen_raw[["gapdh.S"]]

bl_df <- xen_raw[xen_raw$condition == "BL", ]
lb_df <- xen_raw[xen_raw$condition == "LBst", ]

n_bl <- table(factor(bl_df$timepoint, levels = XEN_BL_TP))
n_lb <- table(factor(lb_df$timepoint, levels = XEN_LB_TP))

p_bl_actb <- make_tp_violin(bl_df, "ACTB", "timepoint", "ACTB",
  title    = "Xenopus froglet BL (Lin 2021)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_bl),"=",n_bl), collapse=", "),
                    "  |  actb.L only (actb.S not detected)"),
  tp_levels = XEN_BL_TP, tp_colors = XEN_BL_COL)

p_bl_gapdh <- make_tp_violin(bl_df, "GAPDH", "timepoint", "GAPDH (summed .L+.S)",
  title    = "Xenopus froglet BL (Lin 2021)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_bl),"=",n_bl), collapse=", ")),
  tp_levels = XEN_BL_TP, tp_colors = XEN_BL_COL)

p_lb_actb <- make_tp_violin(lb_df, "ACTB", "timepoint", "ACTB",
  title    = "Xenopus tadpole LBst (Lin 2021)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_lb),"=",n_lb), collapse=", ")),
  tp_levels = XEN_LB_TP, tp_colors = XEN_LB_COL)

p_lb_gapdh <- make_tp_violin(lb_df, "GAPDH", "timepoint", "GAPDH (summed .L+.S)",
  title    = "Xenopus tadpole LBst (Lin 2021)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_lb),"=",n_lb), collapse=", ")),
  tp_levels = XEN_LB_TP, tp_colors = XEN_LB_COL)

rm(xen_raw, bl_df, lb_df); gc()

# =============================================================================
# 3. MOUSE
# =============================================================================

message("Loading Mouse...")
mus <- readRDS(path.expand(MOUSE_RDS))

MUS_TP  <- c("0dpa","11dpa","12dpa","14dpa","17dpa")
MUS_COL <- c("0dpa"="#1B5E20","11dpa"="#2E7D32","12dpa"="#388E3C",
              "14dpa"="#66BB6A","17dpa"="#A5D6A7")

mus_raw <- fetch_genes_split(mus, c("Actb","Gapdh"), meta_cols = "timepoint")
rm(mus); gc()

n_mus <- table(factor(mus_raw$timepoint, levels = MUS_TP))

p_mus_actb <- make_tp_violin(mus_raw, "Actb", "timepoint", "ACTB",
  title    = "Mouse digit tip (Johnson 2020)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_mus),"=",n_mus), collapse=", ")),
  tp_levels = MUS_TP, tp_colors = MUS_COL)

p_mus_gapdh <- make_tp_violin(mus_raw, "Gapdh", "timepoint", "GAPDH",
  title    = "Mouse digit tip (Johnson 2020)",
  subtitle = paste0("n cells: ", paste(paste0(names(n_mus),"=",n_mus), collapse=", ")),
  tp_levels = MUS_TP, tp_colors = MUS_COL)

rm(mus_raw); gc()

# =============================================================================
# 4. COMBINE AND SAVE
# =============================================================================

actg1_page <- (p_axo_actg1 / p_bl_actb / p_lb_actb / p_mus_actb) +
  plot_annotation(
    title    = "Actin expression by condition/timepoint - normalization QC",
    subtitle = "Axolotl: ACTG1 (ACTB absent from annotation). Xenopus/Mouse: ACTB.",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

gapdh_page <- (p_axo_gapdh / p_bl_gapdh / p_lb_gapdh / p_mus_gapdh) +
  plot_annotation(
    title    = "GAPDH expression by condition/timepoint - normalization QC",
    subtitle = "Stable housekeeping expression validates log-normalization across conditions",
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9, color = "grey40"))
  )

out_file <- file.path(RESULTS_DIR, "housekeeping_violin_v2.pdf")
pdf(out_file, width = 10, height = 18)
print(actg1_page)
print(gapdh_page)
dev.off()
message("Saved ", out_file)
