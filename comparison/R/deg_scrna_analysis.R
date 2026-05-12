# =============================================================================
# deg_scrna_analysis.R
#
# Two analyses using BMS vs EtOH DEG lists against the scRNA-seq pseudobulk:
#
#   1. DEG-subset Spearman correlation
#      Same Spearman correlation as the main pipeline, but restricted to genes
#      that are significantly differentially expressed in each bulk condition.
#      Answers: "looking only at genes BMS actually changes, which scRNA-seq
#      state do they most resemble?"
#
#   2. DEG expression landscape heatmap
#      For each DEG gene, plot its mean expression across every scRNA-seq
#      timepoint (z-scored). Rows split by direction (BMS-up vs BMS-down).
#      Answers: "which scRNA-seq states are actively expressing the BMS program?"
#
# Requires bulk_vs_scrna_correlation.R to have been run first so that
# results/pseudobulk_cache.rds exists.
#
# Run from comparison/ directory: source("R/deg_scrna_analysis.R")
#
# Outputs:
#   results/deg_subset_correlation_heatmap.pdf
#   results/deg_expression_landscape.pdf
#   results/deg_subset_correlation_matrix.csv
# =============================================================================

library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

RESULTS_DIR <- "results"
CACHE_FILE  <- file.path(RESULTS_DIR, "pseudobulk_cache.rds")

if (!file.exists(CACHE_FILE))
  stop("Run bulk_vs_scrna_correlation.R first to generate results/pseudobulk_cache.rds")

message("Loading pseudobulk cache...")
cache      <- readRDS(CACHE_FILE)
pseudos    <- cache$pseudos
bulk_mat   <- cache$bulk_mat
rm(cache); gc()

BULK_ORDER <- c("BMS_acute", "BMS_24hrs", "Pretreatment", "Reamputation",
                "EtOH_acute", "EtOH_24hrs", "EtOH_pretreat", "EtOH_reamp")

# =============================================================================
# 1. LOAD DEG TABLES
#    Files are EtOH-vs-BMS comparisons: positive log2FC = higher in EtOH
#    (= lower in BMS = BMS-repressed). We flip sign so positive = BMS-induced.
# =============================================================================

DEG_DIR_2025 <- "~/Desktop/RNA seq BMS Amputation Experiment 2025/deseq2"
DEG_DIR_24HR <- "~/Desktop/24hrs BMS RNA Seq/EtOH-vs-BMS"

deg_paths <- list(
  BMS_acute    = file.path(DEG_DIR_2025, "EtOH_acute-vs-BMS_acute/Significant-DEGs.csv"),
  BMS_24hrs    = file.path(DEG_DIR_24HR, "Significant-DEGs.csv"),
  Pretreatment = file.path(DEG_DIR_2025, "EtOH_pretreat-vs-BMS_pretreat/Significant-DEGs.csv"),
  Reamputation = file.path(DEG_DIR_2025, "EtOH_reamp-vs-BMS_reamp/Significant-DEGs.csv")
)

deg_tables <- lapply(names(deg_paths), function(nm) {
  p  <- path.expand(deg_paths[[nm]])
  df <- read.csv(p, row.names = 1)
  # Flip: positive now means BMS-induced (higher in BMS than EtOH)
  df$log2FC_BMS <- -df$log2FoldChange
  df$direction  <- ifelse(df$log2FC_BMS > 0, "BMS_up", "BMS_down")
  df$condition  <- nm
  df
})
names(deg_tables) <- names(deg_paths)

message("DEG counts per condition:")
for (nm in names(deg_tables))
  message("  ", nm, ": ", nrow(deg_tables[[nm]]),
          " (up: ", sum(deg_tables[[nm]]$direction == "BMS_up"),
          ", down: ", sum(deg_tables[[nm]]$direction == "BMS_down"), ")")

# Union of all DEG gene symbols
all_deg_genes <- unique(unlist(lapply(deg_tables, rownames)))
message("\nTotal unique DEG genes across all conditions: ", length(all_deg_genes))

# Per-gene summary: dominant direction across conditions
deg_summary <- do.call(rbind, lapply(deg_tables, function(df) {
  data.frame(gene = rownames(df), direction = df$direction,
             log2FC_BMS = df$log2FC_BMS, condition = df$condition,
             stringsAsFactors = FALSE)
}))
# For genes in multiple conditions, take majority direction
gene_dir <- tapply(deg_summary$direction, deg_summary$gene, function(x) {
  if (sum(x == "BMS_up") >= sum(x == "BMS_down")) "BMS_up" else "BMS_down"
})
gene_meanfc <- tapply(deg_summary$log2FC_BMS, deg_summary$gene, mean)
gene_conditions <- tapply(deg_summary$condition, deg_summary$gene,
                          function(x) paste(sort(unique(x)), collapse = "+"))

# =============================================================================
# 2. COMBINE ALL scRNA-seq PSEUDOBULK COLUMNS INTO ONE MATRIX
# =============================================================================

all_pb_genes <- Reduce(union, lapply(pseudos, rownames))
all_tp_cols  <- unlist(lapply(pseudos, colnames))

pb_combined <- matrix(NA_real_,
                      nrow     = length(all_pb_genes),
                      ncol     = length(all_tp_cols),
                      dimnames = list(all_pb_genes, all_tp_cols))
for (nm in names(pseudos)) {
  pb <- pseudos[[nm]]
  pb_combined[rownames(pb), colnames(pb)] <- as.matrix(pb)
}

# =============================================================================
# 3. DEG-SUBSET SPEARMAN CORRELATION
# =============================================================================

.spearman_cor_mat <- function(bulk_sub, scrna_sub) {
  cor_out <- matrix(NA_real_, nrow = ncol(bulk_sub), ncol = ncol(scrna_sub),
                    dimnames = list(colnames(bulk_sub), colnames(scrna_sub)))
  for (i in colnames(bulk_sub))
    for (j in colnames(scrna_sub))
      cor_out[i, j] <- cor(bulk_sub[, i], scrna_sub[, j],
                           method = "spearman", use = "complete.obs")
  cor_out
}

message("\n=== Analysis 1: DEG-subset Spearman correlation ===")

# A. Union of all DEGs across conditions
deg_cor_list <- lapply(names(pseudos), function(nm) {
  pb     <- pseudos[[nm]]
  avail  <- intersect(all_deg_genes, intersect(rownames(pb), rownames(bulk_mat)))
  message("  [union] ", nm, ": ", length(avail), " DEGs available")
  if (length(avail) < 5) return(NULL)
  .spearman_cor_mat(bulk_mat[avail, , drop = FALSE],
                    as.matrix(pb)[avail, , drop = FALSE])
})
names(deg_cor_list) <- names(pseudos)
deg_cor_union <- do.call(cbind, Filter(Negate(is.null), deg_cor_list))

message("\nDEG-subset (union) correlation matrix:")
print(round(deg_cor_union, 3))
write.csv(round(deg_cor_union, 4),
          file.path(RESULTS_DIR, "deg_subset_correlation_matrix.csv"))

# B. Per-condition DEG subset (each condition's own DEGs only)
message("\nPer-condition DEG-subset correlations:")
deg_cor_percond <- lapply(names(deg_tables), function(cond) {
  genes <- rownames(deg_tables[[cond]])
  cor_rows <- lapply(names(pseudos), function(nm) {
    pb    <- pseudos[[nm]]
    avail <- intersect(genes, intersect(rownames(pb), rownames(bulk_mat)))
    if (length(avail) < 5) return(NULL)
    row_cond <- if (cond %in% rownames(bulk_mat)) cond else NULL
    if (is.null(row_cond)) return(NULL)
    mat <- .spearman_cor_mat(bulk_mat[avail, row_cond, drop = FALSE],
                             as.matrix(pb)[avail, , drop = FALSE])
    mat
  })
  do.call(cbind, Filter(Negate(is.null), cor_rows))
})
names(deg_cor_percond) <- names(deg_tables)

for (cond in names(deg_cor_percond)) {
  mat <- deg_cor_percond[[cond]]
  if (is.null(mat)) next
  best_tp  <- colnames(mat)[which.max(mat)]
  best_rho <- max(mat, na.rm = TRUE)
  message(sprintf("  %-15s DEGs -> best match: %-30s (rho = %.3f)",
                  cond, best_tp, best_rho))
}

# =============================================================================
# 4. DEG EXPRESSION LANDSCAPE HEATMAP
#    Rows = DEG genes (z-scored across scRNA-seq timepoints)
#    Cols = all scRNA-seq timepoints
#    Row annotation = BMS direction + mean log2FC
# =============================================================================

message("\n=== Analysis 2: DEG expression landscape heatmap ===")

# Subset pseudobulk matrix to DEG genes available in at least one dataset
deg_avail <- intersect(all_deg_genes, rownames(pb_combined))
message("DEG genes with scRNA-seq pseudobulk data: ", length(deg_avail))

deg_pb <- pb_combined[deg_avail, , drop = FALSE]

# Drop columns (timepoints) that are all NA
keep_cols <- colSums(!is.na(deg_pb)) > 0
deg_pb    <- deg_pb[, keep_cols, drop = FALSE]

# Z-score each gene across available timepoints (row-wise)
deg_pb_z <- t(scale(t(deg_pb)))
# Clip extreme z-scores for visual clarity
deg_pb_z <- pmin(pmax(deg_pb_z, -3), 3)

# Row annotations
row_ann <- data.frame(
  Direction = gene_dir[deg_avail],
  Mean_log2FC_BMS = round(gene_meanfc[deg_avail], 2),
  row.names = deg_avail
)

# Column annotations: species label
prefix_to_species <- function(nms) {
  pfx <- sub("_.*", "", nms)
  lut <- c(
    Leigh        = "Axolotl (Leigh)",
    AxoCollab    = "Axolotl (collab)",
    LinBL        = "Xen froglet",
    LinLBst      = "Xen tadpole",
    Mouse        = "Mouse digit tip",
    MouseRegen   = "Mouse regen",
    MouseNonRegen = "Mouse non-regen"
  )
  ifelse(pfx %in% names(lut), lut[pfx], pfx)
}
col_ann <- data.frame(
  Species = prefix_to_species(colnames(deg_pb_z)),
  row.names = colnames(deg_pb_z)
)

# Annotation colors
dir_colors  <- c(BMS_up = "#B2182B", BMS_down = "#2166AC")
n_species   <- length(unique(col_ann$Species))
sp_colors   <- setNames(
  colorRampPalette(brewer.pal(min(n_species, 8), "Set2"))(n_species),
  unique(col_ann$Species)
)
ann_colors  <- list(Direction = dir_colors, Species = sp_colors)

# Sort rows: BMS_up first, then BMS_down, within each group by mean log2FC
row_order <- order(
  row_ann$Direction == "BMS_down",    # BMS_up first (FALSE < TRUE)
  -abs(row_ann$Mean_log2FC_BMS)       # larger FC first
)
deg_pb_z_sorted <- deg_pb_z[row_order, , drop = FALSE]
row_ann_sorted  <- row_ann[row_order, , drop = FALSE]

out_landscape <- file.path(RESULTS_DIR, "deg_expression_landscape.pdf")
n_genes_plot  <- nrow(deg_pb_z_sorted)
pdf(out_landscape,
    width  = max(12, ncol(deg_pb_z_sorted) * 0.55 + 4),
    height = max(8,  n_genes_plot * 0.22 + 3))
pheatmap(
  deg_pb_z_sorted,
  annotation_row    = row_ann_sorted[, "Direction", drop = FALSE],
  annotation_col    = col_ann,
  annotation_colors = ann_colors,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  fontsize_row      = max(5, min(9, 200 / n_genes_plot)),
  fontsize_col      = 8,
  color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks            = seq(-3, 3, length.out = 101),
  border_color      = NA,
  main              = "DEG expression across scRNA-seq pseudobulk timepoints (z-score)",
  angle_col         = 45
)
dev.off()
message("Saved ", out_landscape)

# =============================================================================
# 5. DEG-SUBSET HEATMAP (same layout as main pipeline heatmap)
# =============================================================================

deg_cor_df <- as.data.frame(deg_cor_union) |>
  tibble::rownames_to_column("bulk_condition") |>
  pivot_longer(-bulk_condition, names_to = "scrna_timepoint", values_to = "rho")

deg_cor_df$bulk_condition  <- factor(deg_cor_df$bulk_condition,  levels = BULK_ORDER)
deg_cor_df$scrna_timepoint <- factor(deg_cor_df$scrna_timepoint,
                                     levels = colnames(deg_cor_union))
deg_cor_df$species <- prefix_to_species(as.character(deg_cor_df$scrna_timepoint))
deg_cor_df$species <- factor(deg_cor_df$species,
                             levels = unique(prefix_to_species(colnames(deg_cor_union))))
deg_cor_df$treatment_group <- ifelse(grepl("^EtOH", as.character(deg_cor_df$bulk_condition)),
                                     "EtOH (vehicle)", "BMS-345541")
deg_cor_df$treatment_group <- factor(deg_cor_df$treatment_group,
                                     levels = c("BMS-345541", "EtOH (vehicle)"))

p_deg <- ggplot(deg_cor_df,
                aes(x = scrna_timepoint, y = bulk_condition, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 2.8, color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = median(deg_cor_df$rho, na.rm = TRUE),
    name = "Spearman\nrho"
  ) +
  facet_grid(treatment_group ~ species, scales = "free", space = "free") +
  labs(
    title    = sprintf("DEG-subset Spearman rho (n = %d unique DEG genes)", length(all_deg_genes)),
    subtitle = "Only genes significantly differentially expressed (BMS vs EtOH) used for correlation",
    x = "scRNA-seq timepoint", y = "Bulk condition"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    strip.text       = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 11),
    plot.subtitle    = element_text(size = 8, color = "grey40")
  )

pdf(file.path(RESULTS_DIR, "deg_subset_correlation_heatmap.pdf"),
    width = max(10, ncol(deg_cor_union) * 1.1 + 3),
    height = max(5, nrow(deg_cor_union) * 0.6 + 2))
print(p_deg)
dev.off()
message("Saved deg_subset_correlation_heatmap.pdf")
message("\nDone. Outputs in results/:")
message("  deg_subset_correlation_matrix.csv")
message("  deg_subset_correlation_heatmap.pdf")
message("  deg_expression_landscape.pdf")
