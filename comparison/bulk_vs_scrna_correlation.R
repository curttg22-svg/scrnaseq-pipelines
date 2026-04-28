# =============================================================================
# Bulk RNA-seq vs scRNA-seq pseudobulk correlation
#
# Inputs:
#   Bulk (axolotl whole blastema tissue, normalized counts):
#     ~/Desktop/All RNA normlized count csv files/BMS_normalized_counts.csv
#     ~/Desktop/All RNA normlized count csv files/24 hrs_normalized_counts.csv
#     ~/Desktop/All RNA normlized count csv files/Pretreatment_normalized_counts.csv
#     ~/Desktop/All RNA normlized count csv files/Reamputation_normalized_counts.csv
#
#   scRNA-seq (pseudobulk mean log-normalized expression per timepoint):
#     Leigh 2018 axolotl: ~/Desktop/Axolotl_scrnaseq_project_V3/results/axo_integrated.rds
#     Xenopus (Lin 2021): ~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds
#
# Output: results/bulk_scrna_correlation_heatmap.pdf
#         results/bulk_scrna_correlation_matrix.csv
#
# Method:
#   Spearman correlation between each bulk condition (mean across replicates)
#   and each scRNA-seq pseudobulk timepoint, using shared named gene symbols.
#   LOC IDs (uncharacterized genes) are excluded — they cannot be matched
#   across datasets. Xenopus homeologs (.L/.S) are summed per gene to match
#   what bulk sequencing captures from both subgenome copies.
#
# Interpretation:
#   Higher correlation = the scRNA-seq timepoint's transcriptome most closely
#   resembles your bulk treatment condition at the whole-tissue level.
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

bulk_dir    <- "~/Desktop/All RNA normlized count csv files"
axo_rds     <- "~/Desktop/Axolotl_scrnaseq_project_V3/results/axo_integrated.rds"
xen_rds     <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

WH_REPS <- c("WH_N4","WH_N5","WH_N6")

# =============================================================================
# 1. LOAD BULK DATA
# =============================================================================

message("Loading bulk RNA-seq data...")

load_bulk <- function(path, condition) {
  df <- read.csv(path, row.names = 1, check.names = FALSE)
  # Keep only named genes — LOC IDs cannot be matched across datasets
  named <- df[!grepl("^LOC", rownames(df)), , drop = FALSE]
  # Mean across replicates -> single vector per condition
  setNames(data.frame(rowMeans(named)), condition)
}

bulk_acute <- load_bulk(
  file.path(bulk_dir, "BMS_normalized_counts.csv"),   "BMS_acute")
bulk_24hrs <- load_bulk(
  file.path(bulk_dir, "24 hrs_normalized_counts.csv"), "BMS_24hrs")
bulk_pre   <- load_bulk(
  file.path(bulk_dir, "Pretreatment_normalized_counts.csv"), "Pretreatment")
bulk_reamp <- load_bulk(
  file.path(bulk_dir, "Reamputation_normalized_counts.csv"), "Reamputation")

# Align on shared genes across all bulk conditions
bulk_genes  <- Reduce(intersect, list(rownames(bulk_acute), rownames(bulk_24hrs),
                                      rownames(bulk_pre),   rownames(bulk_reamp)))
bulk_mat <- cbind(
  bulk_acute[bulk_genes, , drop=FALSE],
  bulk_24hrs[bulk_genes, , drop=FALSE],
  bulk_pre[bulk_genes,   , drop=FALSE],
  bulk_reamp[bulk_genes, , drop=FALSE]
)
message("Bulk: ", nrow(bulk_mat), " named genes across 4 conditions")

# =============================================================================
# 2. LEIGH 2018 AXOLOTL — pseudobulk per timepoint
# =============================================================================

message("Computing axolotl pseudobulk...")
axo <- readRDS(axo_rds)
axo <- JoinLayers(axo)

# Pool WH replicates
axo$tp_pooled <- as.character(axo$timepoint)
axo$tp_pooled[axo$tp_pooled %in% WH_REPS] <- "WH"
AXO_TP_ORDER <- c("Intact","WH","EB","MB")
axo$tp_pooled <- factor(axo$tp_pooled, levels = AXO_TP_ORDER)

# Named genes only (filter out Trinity IDs — they start with 'c' + digits)
axo_mat  <- GetAssayData(axo, layer = "data")
named_axo <- !grepl("^c[0-9]", rownames(axo_mat))
axo_mat   <- axo_mat[named_axo, ]

message("Axolotl: ", nrow(axo_mat), " named genes in scRNA-seq object")

axo_pseudo <- do.call(cbind, lapply(AXO_TP_ORDER, function(tp) {
  idx <- which(axo$tp_pooled == tp)
  if (length(idx) == 0) return(NULL)
  setNames(data.frame(Matrix::rowMeans(axo_mat[, idx])),
           paste0("Axo_", tp))
}))

# =============================================================================
# 3. XENOPUS — pseudobulk per timepoint, homeologs summed
# =============================================================================

message("Computing Xenopus pseudobulk...")
xen <- readRDS(xen_rds)
xen <- JoinLayers(xen)

XEN_TP_ORDER <- c("0dpa","3dpa","7-14dpa","14dpa","14-52dpa")

xen_mat  <- GetAssayData(xen, layer = "data")

# Sum .L + .S homeologs per base gene using a sparse indicator matrix multiply.
# This replaces a row-by-row for-loop (2+ hrs) with one matrix operation (<5s).
# Result is mathematically identical: each row of xen_summed is the column sum
# of all homeolog rows for that base gene symbol.
xen_genes    <- rownames(xen_mat)
base_symbols <- toupper(sub("\\.(L|S)(\\.[0-9]+)?$", "", xen_genes))
unique_bases <- unique(base_symbols)
message("Xenopus: summing homeologs for ", length(unique_bases), " base genes")

# Sparse indicator matrix: rows = original genes, cols = base symbols
# Entry [i,j] = 1 if gene i belongs to base symbol j
indicator <- Matrix::sparseMatrix(
  i = seq_along(xen_genes),
  j = match(base_symbols, unique_bases),
  x = 1,
  dims = c(length(xen_genes), length(unique_bases)),
  dimnames = list(xen_genes, unique_bases)
)
# t(indicator) %*% xen_mat sums all homeolog rows into one base-gene row
xen_summed <- Matrix::t(indicator) %*% xen_mat

xen_pseudo <- do.call(cbind, lapply(XEN_TP_ORDER, function(tp) {
  idx <- which(xen$timepoint == tp)
  if (length(idx) == 0) return(NULL)
  setNames(data.frame(Matrix::rowMeans(xen_summed[, idx])),
           paste0("Xen_", tp))
}))

# =============================================================================
# 4. FIND COMMON GENES & COMPUTE SPEARMAN CORRELATIONS
# =============================================================================

message("Finding common genes...")
axo_common <- intersect(rownames(bulk_mat), rownames(axo_pseudo))
xen_common <- intersect(rownames(bulk_mat), rownames(xen_pseudo))
message("Bulk vs Axolotl scRNA-seq: ", length(axo_common), " shared genes")
message("Bulk vs Xenopus scRNA-seq: ", length(xen_common), " shared genes")

# Spearman correlation function
spearman_cor <- function(bulk_mat, scrna_pseudo, common_genes) {
  b <- bulk_mat[common_genes, , drop=FALSE]
  s <- scrna_pseudo[common_genes, , drop=FALSE]
  cor_mat <- matrix(NA, ncol(b), ncol(s),
                    dimnames = list(colnames(b), colnames(s)))
  for (i in colnames(b))
    for (j in colnames(s))
      cor_mat[i,j] <- cor(b[,i], s[,j], method="spearman", use="complete.obs")
  cor_mat
}

cor_axo <- spearman_cor(bulk_mat, axo_pseudo, axo_common)
cor_xen <- spearman_cor(bulk_mat, xen_pseudo, xen_common)

# Combine into one matrix: axolotl timepoints then Xenopus timepoints
cor_all <- cbind(cor_axo, cor_xen)

message("\nCorrelation summary:")
print(round(cor_all, 3))

# Save CSV
write.csv(round(cor_all, 4),
          file.path(results_dir, "bulk_scrna_correlation_matrix.csv"))
message("Saved bulk_scrna_correlation_matrix.csv")

# =============================================================================
# 5. HEATMAP
# =============================================================================

# Reshape for ggplot
cor_df <- as.data.frame(cor_all) |>
  tibble::rownames_to_column("bulk_condition") |>
  pivot_longer(-bulk_condition, names_to = "scrna_timepoint", values_to = "rho")

# Order axes
BULK_ORDER <- c("BMS_acute","BMS_24hrs","Pretreatment","Reamputation")
SCRNA_ORDER <- c(
  paste0("Axo_", c("Intact","WH","EB","MB")),
  paste0("Xen_", c("0dpa","3dpa","7-14dpa","14dpa","14-52dpa"))
)
cor_df$bulk_condition  <- factor(cor_df$bulk_condition,  levels = BULK_ORDER)
cor_df$scrna_timepoint <- factor(cor_df$scrna_timepoint, levels = SCRNA_ORDER)

# Add species label for x-axis grouping
cor_df$species <- ifelse(grepl("^Axo_", cor_df$scrna_timepoint),
                          "Axolotl (Leigh 2018)", "Xenopus (Lin 2021)")

p <- ggplot(cor_df, aes(x = scrna_timepoint, y = bulk_condition, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(rho, 3)), size = 3, color = "black") +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#B2182B",
    midpoint = median(cor_df$rho, na.rm = TRUE),
    name     = "Spearman\ncorrelation"
  ) +
  facet_grid(. ~ species, scales = "free_x", space = "free_x") +
  labs(
    title    = "Bulk RNA-seq vs scRNA-seq pseudobulk: Spearman correlation",
    subtitle = paste0(
      "Axolotl scRNA-seq: ", length(axo_common), " shared genes  |  ",
      "Xenopus scRNA-seq: ", length(xen_common), " shared genes\n",
      "Bulk: whole axolotl blastema (BMS experiment); scRNA-seq: whole-timepoint pseudobulk"
    ),
    x = "scRNA-seq timepoint",
    y = "Bulk condition"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    legend.position  = "right"
  )

pdf(file.path(results_dir, "bulk_scrna_correlation_heatmap.pdf"),
    width = 12, height = 5)
print(p)
dev.off()
message("Saved bulk_scrna_correlation_heatmap.pdf")

# =============================================================================
# 6. BEST-MATCH SUMMARY
# =============================================================================

message("\n--- Best-matching scRNA-seq timepoint per bulk condition ---")
for (cond in BULK_ORDER) {
  if (!cond %in% rownames(cor_all)) next
  best <- names(which.max(cor_all[cond, ]))
  message(sprintf("  %-15s -> %s  (rho = %.3f)", cond, best, cor_all[cond, best]))
}
