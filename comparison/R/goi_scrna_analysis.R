# =============================================================================
# goi_scrna_analysis.R
#
# Genes-of-interest (GOI) analysis against the scRNA-seq pseudobulk compendium.
#
# Two gene sets are analysed independently:
#
#   GOI_ALL  — all 16 genes of interest
#   GOI_HH   — Hedgehog signaling components (ligands, receptor, transducer,
#               transcription factors, and canonical target/feedback genes)
#
# For each gene set:
#   1. Spearman correlation using only those genes (shows which scRNA-seq state
#      has the most similar expression profile for these specific genes)
#   2. Expression dot/heatmap: mean log-norm expression of each gene across
#      every scRNA-seq timepoint — raw values, not z-scored, so absolute
#      expression level is visible
#
# Requires bulk_vs_scrna_correlation.R to have been run first so that
# results/pseudobulk_cache.rds exists.
#
# Run from comparison/ directory: source("R/goi_scrna_analysis.R")
#
# Outputs:
#   results/goi_all_correlation_heatmap.pdf
#   results/goi_all_expression_heatmap.pdf
#   results/goi_hh_correlation_heatmap.pdf
#   results/goi_hh_expression_heatmap.pdf
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
cache    <- readRDS(CACHE_FILE)
pseudos  <- cache$pseudos
bulk_mat <- cache$bulk_mat
rm(cache); gc()

BULK_ORDER <- c("BMS_acute", "BMS_24hrs", "Pretreatment", "Reamputation",
                "EtOH_acute", "EtOH_24hrs", "EtOH_pretreat", "EtOH_reamp")

# =============================================================================
# GENE SETS
# =============================================================================

GOI_ALL <- c(
  "C3", "LDLR", "MSMO1", "CYP51A1", "FGFBP1", "FDPS",
  "GDF5", "IL11", "MSX1", "CLDN5", "CYTL1", "FREM1",
  "GLI1", "GREM1", "PTCH1", "PTCH2"
)

# Core Hedgehog pathway: ligands, receptor, transducer, GLI TFs, and
# canonical transcriptional targets / feedback regulators
GOI_HH <- c(
  "SHH", "IHH", "DHH",          # ligands
  "PTCH1", "PTCH2",             # receptors (PTCH1 is also a transcriptional target)
  "SMO",                         # signal transducer
  "GLI1", "GLI2", "GLI3",       # transcription factors
  "HHIP",                        # pathway antagonist / HH target
  "BOC", "GAS1", "CDON"         # co-receptors
)

# =============================================================================
# HELPERS
# =============================================================================

prefix_to_species <- function(nms) {
  pfx <- sub("_.*", "", nms)
  lut <- c(
    Leigh         = "Axolotl (Leigh)",
    AxoCollab     = "Axolotl (collab)",
    LinBL         = "Xen froglet",
    LinLBst       = "Xen tadpole",
    Mouse         = "Mouse digit tip",
    MouseRegen    = "Mouse regen",
    MouseNonRegen = "Mouse non-regen"
  )
  ifelse(pfx %in% names(lut), lut[pfx], pfx)
}

.spearman_cor_mat <- function(bulk_sub, scrna_sub) {
  cor_out <- matrix(NA_real_, nrow = ncol(bulk_sub), ncol = ncol(scrna_sub),
                    dimnames = list(colnames(bulk_sub), colnames(scrna_sub)))
  for (i in colnames(bulk_sub))
    for (j in colnames(scrna_sub))
      cor_out[i, j] <- cor(bulk_sub[, i], scrna_sub[, j],
                           method = "spearman", use = "complete.obs")
  cor_out
}

# Build combined scRNA-seq pseudobulk matrix (all datasets, all timepoints)
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
# MAIN ANALYSIS FUNCTION — runs both analyses for one gene set
# =============================================================================

run_goi_analysis <- function(gene_set, label, label_short) {

  message("\n", strrep("=", 60))
  message("Gene set: ", label, " (", length(gene_set), " genes)")
  message(strrep("=", 60))

  # ---- 1. Which genes are actually present in each dataset? -----------------
  avail_per_dataset <- lapply(names(pseudos), function(nm) {
    intersect(gene_set, rownames(pseudos[[nm]]))
  })
  names(avail_per_dataset) <- names(pseudos)

  avail_in_bulk    <- intersect(gene_set, rownames(bulk_mat))
  avail_in_scrna   <- intersect(gene_set, rownames(pb_combined))
  missing_from_all <- setdiff(gene_set, union(avail_in_bulk, avail_in_scrna))

  message("Genes present in bulk:  ", paste(avail_in_bulk, collapse=", "))
  message("Genes present in scRNA: ", paste(avail_in_scrna, collapse=", "))
  if (length(missing_from_all) > 0)
    message("Not found anywhere:     ", paste(missing_from_all, collapse=", "))

  # ---- 2. Spearman correlation (GOI subset) ---------------------------------
  cor_list <- lapply(names(pseudos), function(nm) {
    pb    <- pseudos[[nm]]
    avail <- intersect(gene_set, intersect(rownames(pb), rownames(bulk_mat)))
    message("  ", nm, ": ", length(avail), " genes used for correlation")
    if (length(avail) < 3) return(NULL)
    .spearman_cor_mat(bulk_mat[avail, , drop = FALSE],
                      as.matrix(pb)[avail, , drop = FALSE])
  })
  names(cor_list) <- names(pseudos)
  cor_mat <- do.call(cbind, Filter(Negate(is.null), cor_list))

  message("\nCorrelation matrix (", label_short, "):")
  print(round(cor_mat, 3))

  message("\nBest-matching timepoint per bulk condition:")
  for (cond in BULK_ORDER) {
    if (!cond %in% rownames(cor_mat)) next
    row  <- cor_mat[cond, ]
    best <- names(which.max(row))
    message(sprintf("  %-15s -> %-30s (rho = %.3f)", cond, best, row[best]))
  }

  write.csv(round(cor_mat, 4),
            file.path(RESULTS_DIR, paste0("goi_", label_short, "_correlation_matrix.csv")))

  # Correlation heatmap
  cor_df <- as.data.frame(cor_mat) |>
    tibble::rownames_to_column("bulk_condition") |>
    pivot_longer(-bulk_condition, names_to = "scrna_timepoint", values_to = "rho")

  cor_df$bulk_condition  <- factor(cor_df$bulk_condition,  levels = BULK_ORDER)
  cor_df$scrna_timepoint <- factor(cor_df$scrna_timepoint, levels = colnames(cor_mat))
  cor_df$species         <- prefix_to_species(as.character(cor_df$scrna_timepoint))
  cor_df$species         <- factor(cor_df$species,
                                   levels = unique(prefix_to_species(colnames(cor_mat))))
  cor_df$treatment_group <- ifelse(grepl("^EtOH", as.character(cor_df$bulk_condition)),
                                   "EtOH (vehicle)", "BMS-345541")
  cor_df$treatment_group <- factor(cor_df$treatment_group,
                                   levels = c("BMS-345541", "EtOH (vehicle)"))

  p_cor <- ggplot(cor_df, aes(x = scrna_timepoint, y = bulk_condition, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", rho)), size = 2.8, color = "black") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = median(cor_df$rho, na.rm = TRUE),
      name = "Spearman\nrho"
    ) +
    facet_grid(treatment_group ~ species, scales = "free", space = "free") +
    labs(
      title    = paste0(label, ": Spearman rho vs scRNA-seq pseudobulk"),
      subtitle = paste0("n = ", length(avail_in_scrna), " genes found in scRNA-seq  |  ",
                        length(missing_from_all), " genes not detected in any dataset"),
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

  cor_pdf <- file.path(RESULTS_DIR, paste0("goi_", label_short, "_correlation_heatmap.pdf"))
  pdf(cor_pdf, width = max(10, ncol(cor_mat) * 1.1 + 3),
      height = max(5, nrow(cor_mat) * 0.6 + 2))
  print(p_cor)
  dev.off()
  message("Saved ", cor_pdf)

  # ---- 3. Expression heatmap across scRNA-seq timepoints --------------------
  # Raw log-norm expression (not z-scored) so absolute level is visible.
  # Rows = GOI genes, columns = scRNA-seq timepoints.
  goi_avail <- intersect(gene_set, rownames(pb_combined))
  goi_pb    <- pb_combined[goi_avail, , drop = FALSE]

  # Drop all-NA columns
  keep_cols <- colSums(!is.na(goi_pb)) > 0
  goi_pb    <- goi_pb[, keep_cols, drop = FALSE]

  # Column annotation
  col_ann <- data.frame(
    Species = prefix_to_species(colnames(goi_pb)),
    row.names = colnames(goi_pb)
  )
  n_sp      <- length(unique(col_ann$Species))
  sp_colors <- setNames(
    colorRampPalette(brewer.pal(min(n_sp, 8), "Set2"))(n_sp),
    unique(col_ann$Species)
  )

  expr_pdf <- file.path(RESULTS_DIR, paste0("goi_", label_short, "_expression_heatmap.pdf"))
  n_genes  <- nrow(goi_pb)
  pdf(expr_pdf,
      width  = max(12, ncol(goi_pb) * 0.55 + 4),
      height = max(5,  n_genes * 0.45 + 3))
  pheatmap(
    goi_pb,
    annotation_col    = col_ann,
    annotation_colors = list(Species = sp_colors),
    cluster_rows      = n_genes > 2,   # cluster if >2 genes
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 10,
    fontsize_col      = 8,
    color             = colorRampPalette(c("#F7F7F7", "#FEE08B", "#D73027"))(100),
    border_color      = "grey80",
    main              = paste0(label, ": mean log-norm expression across scRNA-seq timepoints"),
    angle_col         = 45,
    na_col            = "grey70"
  )
  dev.off()
  message("Saved ", expr_pdf)
}

# =============================================================================
# RUN BOTH GENE SETS
# =============================================================================

run_goi_analysis(GOI_ALL, label = "Genes of interest (all 16)", label_short = "all")
run_goi_analysis(GOI_HH,  label = "Hedgehog signaling components", label_short = "hh")

message("\nDone. Outputs in results/:")
message("  goi_all_correlation_matrix.csv")
message("  goi_all_correlation_heatmap.pdf")
message("  goi_all_expression_heatmap.pdf")
message("  goi_hh_correlation_matrix.csv")
message("  goi_hh_correlation_heatmap.pdf")
message("  goi_hh_expression_heatmap.pdf")
