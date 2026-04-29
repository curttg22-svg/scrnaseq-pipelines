# =============================================================================
# plot_correlation_scatters.R
#
# For each bulk condition, generates a scatter plot of bulk expression vs
# the best-matching scRNA-seq pseudobulk timepoint.
#
# Run from comparison/ directory: source("R/plot_correlation_scatters.R")
#
# Output: results/correlation_scatter_bestmatch.pdf
# =============================================================================

library(Seurat)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(patchwork)

source("R/utils_harmonize.R")

BULK_DIR  <- "~/Desktop/All RNA normlized count csv files"
AXO_RDS   <- "~/Desktop/Axolotl_scrnaseq_project_V3/results/axo_integrated.rds"
XEN_RDS   <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
MOUSE_RDS <- "~/Desktop/scrnaseq-pipelines/mouse-digit/results/mouse_merged.rds"
RESULTS_DIR <- "results"

WH_REPS    <- c("WH_N4","WH_N5","WH_N6")
BULK_ORDER <- c("BMS_acute","BMS_24hrs","Pretreatment","Reamputation")

# Genes to label in scatter plots — canonical limb/blastema markers
LABEL_GENES <- c(
  "COL1A1","COL1A2","COL3A1","VIM","DCN","FN1",    # CT/fibroblast
  "COL2A1","SOX9","ACAN",                            # cartilage
  "KRT8","KRT18","EPCAM","KRT5",                     # epithelial
  "PTPRC","CD68","S100A8",                           # immune
  "PECAM1","CDH5","VWF",                             # endothelial
  "ACTA2","TAGLN","MYH11",                           # smooth muscle
  "SOX10","S100B",                                   # neural/Schwann
  "MKI67","TOP2A",                                   # proliferating
  "HBA1","HBB"                                       # erythrocyte
)

# =============================================================================
# 1. HELPERS (duplicated from main script for standalone use)
# =============================================================================

.slim_for_pseudobulk <- function(obj) {
  for (l in Layers(obj)[startsWith(Layers(obj), "counts.")])
    LayerData(obj, layer = l) <- NULL
  for (l in Layers(obj)[startsWith(Layers(obj), "scale.data")])
    LayerData(obj, layer = l) <- NULL
  gc()
  JoinLayers(obj)
}

.mean_per_tp <- function(expr_mat, tp_factor, tp_order, prefix) {
  parts <- lapply(tp_order, function(tp) {
    idx <- which(tp_factor == tp)
    if (length(idx) == 0) return(NULL)
    m <- Matrix::rowMeans(expr_mat[, idx, drop = FALSE])
    data.frame(m, row.names = names(m))
  })
  parts <- Filter(Negate(is.null), parts)
  result <- do.call(cbind, parts)
  colnames(result) <- paste0(prefix, "_", tp_order[seq_along(parts)])
  result
}

# =============================================================================
# 2. BULK DATA
# =============================================================================

message("Loading bulk RNA-seq data...")
read_bulk <- function(fname, condition) {
  df <- read.csv(file.path(BULK_DIR, fname), row.names = 1, check.names = FALSE)
  mat <- as.matrix(df)
  genes <- toupper(rownames(mat))
  genes <- genes[!grepl("^LOC[0-9]", genes)]
  mat   <- mat[!grepl("^LOC[0-9]", toupper(rownames(mat))), , drop = FALSE]
  rownames(mat) <- genes
  rowMeans(mat)
}

bulk_vec <- list(
  BMS_acute    = read_bulk("BMS_normalized_counts.csv",      "BMS_acute"),
  BMS_24hrs    = read_bulk("24 hrs_normalized_counts.csv",   "BMS_24hrs"),
  Pretreatment = read_bulk("Pretreatment_normalized_counts.csv", "Pretreatment"),
  Reamputation = read_bulk("Reamputation_normalized_counts.csv", "Reamputation")
)

# =============================================================================
# 3. PSEUDOBULK (all three species — recomputed for scatter access)
# =============================================================================

message("Computing pseudobulk: Axolotl (Leigh 2018)...")
axo  <- .slim_for_pseudobulk(readRDS(AXO_RDS))
expr <- GetAssayData(axo, layer = "data")
tp   <- as.character(axo$timepoint); tp[tp %in% WH_REPS] <- "WH"
tp_order <- c("Intact","WH","EB","MB")
pb_axo <- .mean_per_tp(expr, factor(tp, levels = tp_order), tp_order, prefix = "Leigh")
pb_axo <- aggregate_to_hgnc(as.matrix(pb_axo), species = "human",
                              filter_fn = function(g) !grepl("^c[0-9]", g))
rm(axo, expr); gc()

message("Computing pseudobulk: Xenopus (Lin 2021)...")
xen      <- .slim_for_pseudobulk(readRDS(XEN_RDS))
expr     <- GetAssayData(xen, layer = "data")
xen_genes    <- rownames(expr)
base_symbols <- toupper(sub("\\.(L|S)(\\.[0-9]+)?$", "", xen_genes))
unique_bases <- unique(base_symbols)
indicator <- sparseMatrix(
  i = seq_along(xen_genes), j = match(base_symbols, unique_bases), x = 1,
  dims = c(length(xen_genes), length(unique_bases)),
  dimnames = list(xen_genes, unique_bases)
)
expr_summed <- t(indicator) %*% expr
bl_idx <- which(xen$condition == "BL")
lb_idx <- which(xen$condition == "LBst")
pb_bl  <- .mean_per_tp(expr_summed[, bl_idx],
                        factor(xen$timepoint[bl_idx],
                               levels = c("0dpa","3dpa","7-14dpa","14dpa","14-52dpa")),
                        c("0dpa","3dpa","7-14dpa","14dpa","14-52dpa"), prefix = "LinBL")
pb_lb  <- .mean_per_tp(expr_summed[, lb_idx],
                        factor(xen$timepoint[lb_idx], levels = c("NF50","NF51","NF52")),
                        c("NF50","NF51","NF52"), prefix = "LinLBst")
pb_xen <- cbind(pb_bl, pb_lb)
rm(xen, expr, expr_summed, indicator); gc()

message("Computing pseudobulk: Mouse digit tip...")
mus     <- .slim_for_pseudobulk(readRDS(MOUSE_RDS))
expr    <- GetAssayData(mus, layer = "data")
tp_order_mus <- c("0dpa","11dpa","12dpa","14dpa","17dpa")
pb_mus  <- .mean_per_tp(expr, factor(mus$timepoint, levels = tp_order_mus),
                         tp_order_mus, prefix = "Mouse")
pb_mus  <- aggregate_to_hgnc(as.matrix(pb_mus), species = "mouse")
rm(mus, expr); gc()

# Combine all pseudobulk timepoints into one lookup
all_pb <- list(Leigh = pb_axo, Xenopus = as.matrix(pb_xen), Mouse = pb_mus)

# =============================================================================
# 4. FIND BEST MATCH PER BULK CONDITION
# =============================================================================

message("Finding best matches...")

# Build full combined matrix across all scRNA-seq timepoints
all_genes  <- Reduce(intersect, lapply(all_pb, rownames))
pb_combined <- do.call(cbind, lapply(all_pb, function(m) m[all_genes, ]))

best_matches <- lapply(BULK_ORDER, function(cond) {
  b_vec <- bulk_vec[[cond]]
  shared <- intersect(names(b_vec), all_genes)
  b <- b_vec[shared]
  rhos <- apply(pb_combined[shared, ], 2, function(s)
    cor(b, s, method = "spearman", use = "complete.obs"))
  best_col <- names(which.max(rhos))
  list(bulk_cond = cond, best_col = best_col, rho = max(rhos),
       shared = shared, b = b, s = pb_combined[shared, best_col])
})
names(best_matches) <- BULK_ORDER

# =============================================================================
# 5. SCATTER PLOTS
# =============================================================================

make_scatter <- function(bm) {
  df <- data.frame(
    bulk   = bm$b,
    scrna  = bm$s,
    gene   = names(bm$b)
  )
  df$label <- ifelse(df$gene %in% LABEL_GENES, df$gene, NA)

  # Friendly axis labels
  bulk_label  <- bm$bulk_cond
  scrna_label <- bm$best_col

  ggplot(df, aes(x = log1p(bulk), y = scrna)) +
    geom_point(alpha = 0.25, size = 0.6, color = "#546E7A") +
    geom_smooth(method = "lm", se = FALSE, color = "#D32F2F",
                linewidth = 0.8, linetype = "solid") +
    geom_label_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                     fill = "white", label.padding = 0.15,
                     segment.color = "#424242", segment.size = 0.3,
                     na.rm = TRUE) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("rho = %.3f\nn = %d genes", bm$rho, nrow(df)),
             hjust = -0.1, vjust = 1.3, size = 3.5, fontface = "bold") +
    labs(
      title    = sprintf("%s vs %s", bulk_label, scrna_label),
      x        = "Bulk RNA-seq (log1p normalized counts)",
      y        = "scRNA-seq pseudobulk (mean log-norm)"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

plots <- lapply(best_matches, make_scatter)

combined <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
combined <- combined +
  plot_annotation(
    title    = "Bulk RNA-seq vs best-matching scRNA-seq pseudobulk",
    subtitle = "Each panel: one bulk condition vs its best-matching scRNA-seq timepoint (Spearman)",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 10))
  )

out_file <- file.path(RESULTS_DIR, "correlation_scatter_bestmatch.pdf")
pdf(out_file, width = 14, height = 11)
print(combined)
dev.off()
message("Saved ", out_file)

# Print summary
message("\nBest match per bulk condition:")
for (nm in BULK_ORDER) {
  bm <- best_matches[[nm]]
  message(sprintf("  %-15s -> %-20s (rho = %.3f, n = %d genes)",
          bm$bulk_cond, bm$best_col, bm$rho, length(bm$shared)))
}
