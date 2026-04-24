# =============================================================================
# 06 — Pseudobulk DESeq2: WH vs Intact
# Input:  results/axo_integrated.rds
# Output: results/pseudobulk/  — CSV tables + PDF figures
#
# Pseudobulk approach: aggregate raw counts per (cell_type x replicate),
# then run DESeq2 with replicate as the unit of observation.
#
# Why pseudobulk over cell-level Wilcoxon:
#   - Cells from the same replicate are not independent observations
#   - Wilcoxon on all cells gives inflated significance (flat p-value ceiling)
#   - DESeq2 on aggregated counts respects biological replicate structure
#
# WH has 3 replicates (N4/N5/N6) vs Intact (1 replicate = SCP422).
# This is the only comparison in the Leigh 2018 dataset with true replicates.
# EB and MB are single samples — pseudobulk is not valid for those comparisons.
# For EB and MB use 07_exploratory_de.R (cell-level Wilcoxon, treated as
# exploratory only due to absence of biological replicates).
#
# Reference: Squair et al. 2021, Nat Commun — pseudobulk outperforms
#            cell-level methods for differential abundance and expression.
#
# ANALYSIS NOTES:
#   CIRBP (Cold-Inducible RNA Binding Protein) is consistently the top
#   downregulated gene across every cell type in WH vs Intact. Two possible
#   explanations:
#   (1) Real global wound response — CIRBP is known to drop rapidly after
#       tissue injury as cells transition from resting to active state.
#   (2) Batch/processing artifact — SCP422 (Intact) and SCP489 (WH) are
#       separate studies; any difference in sample processing time or
#       temperature would be captured by CIRBP due to its stress sensitivity.
#   With only one Intact sample these two explanations cannot be distinguished.
#   CIRBP results should be noted but not over-interpreted as cell-type-specific
#   biology. Flag in any downstream analysis or publication.
#
#   Cell type harmonization: SCP422 and SCP489 used different cluster naming
#   conventions. CT_HARMONIZE maps equivalent cell types across datasets.
#   See table below for mappings used.
# =============================================================================

library(Seurat)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(patchwork)

results_dir <- "results"
pb_dir      <- file.path(results_dir, "pseudobulk")
dir.create(pb_dir, showWarnings = FALSE, recursive = TRUE)

WH_REPS  <- c("WH_N4", "WH_N5", "WH_N6")

# =============================================================================
# CELL TYPE HARMONIZATION MAP
# SCP422 (Intact) and SCP489 (WH) used different naming conventions.
# This table maps equivalent cell types across the two datasets.
# Only matched pairs are used for pseudobulk DE.
# =============================================================================

CT_HARMONIZE <- data.frame(
  intact_label = c(
    "SSCs",
    "Fibroblast",
    "Tcells",
    "Early_B",
    "RBCs",
    "Myeloid",
    "BasalEpidermis",
    "Intermediate_epidermis",
    "Proliferating_Epidermis",
    "Epidermal_Langerhans"
  ),
  wh_label = c(
    "SSC",
    "Fibroblast-like_blastema_",
    "T_cell",
    "Early_B_cell",
    "Erythrocyte",
    "Recruited_Macrophage",
    "Basal WE",
    "Intermediate_WE_",
    "Intermediate_WE_",
    "Dendritic_Cell"
  ),
  harmonized = c(
    "SSC_fibroblast",
    "Fibroblast_blastema",
    "T_cell",
    "Early_B_cell",
    "Erythrocyte_RBC",
    "Myeloid_macrophage",
    "Basal_epidermis",
    "Intermediate_epidermis",
    "Intermediate_epidermis",
    "Dendritic_Langerhans"
  ),
  stringsAsFactors = FALSE
)

GOI <- c(
  "GLI1","PTCH1","PTCH2","SHH","IHH","DHH",
  "GDF5","MSX1","SALL1","GREM1",
  "FREM1","CYTL1","CLDN5","IL11","FGFBP1",
  "C3","LDLR","ALOX15B","MSMO1","CYP51A1","FDPS"
)

# =============================================================================
# 1. LOAD & SUBSET TO WH + INTACT ONLY
# =============================================================================

message("Loading axo_integrated.rds...")
axo <- readRDS(file.path(results_dir, "axo_integrated.rds"))
axo <- JoinLayers(axo)

# Keep only Intact and WH cells; label WH replicate
axo$replicate <- as.character(axo$timepoint)
axo_wh_intact <- subset(axo, timepoint %in% c("Intact", WH_REPS))

# Assign condition: Intact vs WH
axo_wh_intact$condition <- ifelse(
  axo_wh_intact$timepoint == "Intact", "Intact", "WH"
)

# Use replicate ID as the pseudobulk sample label
# Intact = one sample ("Intact"); WH = three samples (N4, N5, N6)
axo_wh_intact$pb_sample <- axo_wh_intact$replicate

message("Cells: ", ncol(axo_wh_intact))
message("Samples:")
print(table(axo_wh_intact$pb_sample, axo_wh_intact$condition))

# Apply harmonized cell type labels
axo_wh_intact$ct_harmonized <- NA_character_

# Map Intact labels
for (i in seq_len(nrow(CT_HARMONIZE))) {
  idx <- axo_wh_intact$condition == "Intact" &
         axo_wh_intact$paper_cluster == CT_HARMONIZE$intact_label[i]
  axo_wh_intact$ct_harmonized[idx] <- CT_HARMONIZE$harmonized[i]
}
# Map WH labels
for (i in seq_len(nrow(CT_HARMONIZE))) {
  idx <- axo_wh_intact$condition == "WH" &
         axo_wh_intact$paper_cluster == CT_HARMONIZE$wh_label[i]
  axo_wh_intact$ct_harmonized[idx] <- CT_HARMONIZE$harmonized[i]
}

# Drop cells with no harmonized label
axo_wh_intact <- subset(axo_wh_intact,
                         !is.na(ct_harmonized))

message("Cells with harmonized labels: ", ncol(axo_wh_intact))
message("Harmonized cell type counts:")
print(table(axo_wh_intact$ct_harmonized, axo_wh_intact$condition))

# =============================================================================
# 2. AGGREGATE COUNTS PER (CELL_TYPE x SAMPLE) — PSEUDOBULK
# =============================================================================

message("\nAggregating counts (pseudobulk)...")

counts_mat <- GetAssayData(axo_wh_intact, layer = "counts")
ct_vec     <- as.character(axo_wh_intact$ct_harmonized)
samp_vec   <- as.character(axo_wh_intact$pb_sample)
cond_vec   <- as.character(axo_wh_intact$condition)

# Get harmonized cell types with enough cells in both conditions
ct_table <- table(ct_vec, cond_vec)
valid_ct <- rownames(ct_table)[
  ct_table[, "Intact"] >= 10 &
  ct_table[, "WH"]     >= 30
]
message("Cell types with sufficient cells: ", length(valid_ct))
print(valid_ct)

# Aggregate: sum raw counts per (cell_type x sample)
aggregate_counts <- function(ct) {
  samples <- unique(samp_vec)
  pb_list <- lapply(samples, function(s) {
    idx <- which(ct_vec == ct & samp_vec == s)
    if (length(idx) == 0) return(NULL)
    Matrix::rowSums(counts_mat[, idx, drop = FALSE])
  })
  names(pb_list) <- samples
  pb_list <- Filter(Negate(is.null), pb_list)
  mat <- do.call(cbind, pb_list)
  round(mat)
}

# =============================================================================
# 3. RUN DESeq2 PER CELL TYPE
# =============================================================================

message("\nRunning DESeq2 per cell type...")

all_de <- list()

for (ct in valid_ct) {
  message("  Processing: ", ct)

  pb_mat <- aggregate_counts(ct)
  if (ncol(pb_mat) < 3) {
    message("    Skipping — fewer than 3 samples after aggregation")
    next
  }

  # Build sample metadata
  sample_names <- colnames(pb_mat)
  col_data <- data.frame(
    sample    = sample_names,
    condition = ifelse(sample_names == "Intact", "Intact", "WH"),
    row.names = sample_names
  )
  col_data$condition <- factor(col_data$condition, levels = c("Intact", "WH"))

  # Filter lowly expressed genes: keep genes with > 5 counts in at least 2 samples
  # Filter 1: must have > 5 counts in at least 2 samples
  keep1 <- rowSums(pb_mat > 5) >= 2
  # Filter 2: must have > 0 counts in the Intact sample (removes zero-inflation artifacts)
  intact_col <- colnames(pb_mat)[colnames(pb_mat) == "Intact"]
  keep2 <- if (length(intact_col) > 0) pb_mat[, intact_col] > 0 else TRUE
  pb_mat_filt <- pb_mat[keep1 & keep2, ]
  message("    Genes after filter: ", nrow(pb_mat_filt),
          " / ", nrow(pb_mat))

  if (nrow(pb_mat_filt) < 100) {
    message("    Skipping — too few genes after filtering")
    next
  }

  # DESeq2
  dds <- tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = pb_mat_filt,
      colData   = col_data,
      design    = ~ condition
    )
    DESeq(dds, quiet = TRUE)
  }, error = function(e) {
    message("    DESeq2 error: ", e$message)
    NULL
  })

  if (is.null(dds)) next

  res <- results(dds, contrast = c("condition", "WH", "Intact"),
                 alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene      <- rownames(res_df)
  res_df$cell_type <- ct

  all_de[[ct]] <- res_df

  ct_label <- gsub("[^A-Za-z0-9_]", "_", ct)
  write.csv(res_df,
            file.path(pb_dir, paste0(ct_label, "_WH_vs_Intact_pseudobulk.csv")),
            row.names = FALSE)

  n_sig <- sum(res_df$padj < 0.05 & !is.na(res_df$padj))
  message("    Significant DEGs (padj<0.05): ", n_sig)
}

message("\nCompleted ", length(all_de), " cell types")

# =============================================================================
# 4. VOLCANO PLOTS (one per cell type)
# =============================================================================

message("\nGenerating volcano plots...")

make_pb_volcano <- function(res_df, title, goi = GOI) {
  res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
  res_df$sig <- "NS"
  res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "Up"
  res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
  res_df$sig <- factor(res_df$sig, levels = c("Up","Down","NS"))

  # Cap axes
  res_df$log10p  <- pmin(-log10(res_df$padj + 1e-300), 30)
  res_df$fc_plot <- pmax(pmin(res_df$log2FoldChange, 8), -8)

  # Label: only clean gene symbols
  res_named <- res_df |> filter(!grepl("^c[0-9]|-g[0-9]|-i[0-9]|\\^", gene))
  top_up    <- res_named |> filter(sig == "Up")   |> arrange(desc(log2FoldChange)) |> head(8)
  top_down  <- res_named |> filter(sig == "Down")  |> arrange(log2FoldChange)      |> head(8)
  goi_hits  <- res_df    |> filter(gene %in% goi, padj < 0.05)
  label_df  <- bind_rows(top_up, top_down, goi_hits) |> distinct(gene, .keep_all = TRUE)

  n_up   <- sum(res_df$sig == "Up")
  n_down <- sum(res_df$sig == "Down")

  ggplot(res_df, aes(x = fc_plot, y = log10p, color = sig)) +
    geom_point(size = 0.9, alpha = 0.5) +
    geom_point(data = label_df, size = 2, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = label_df, aes(label = gene),
      size = 2.8, max.overlaps = 15,
      box.padding = 0.4, color = "grey20"
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    scale_color_manual(
      values = c("Up" = "#B71C1C", "Down" = "#1565C0", "NS" = "grey75"),
      name = NULL
    ) +
    labs(title    = title,
         subtitle = paste0("Up: ", n_up, "  |  Down: ", n_down,
                           "  |  padj<0.05, |log2FC|>1"),
         x = "log2 fold-change WH vs Intact (capped at +/-8)",
         y = "-log10 adjusted p-value (capped at 30)") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")
}

pdf(file.path(pb_dir, "pseudobulk_volcano_plots.pdf"), width = 7, height = 5)
for (ct in names(all_de)) {
  print(make_pb_volcano(all_de[[ct]],
                        title = paste0(ct, " - WH vs Intact (pseudobulk DESeq2)")))
}
dev.off()
message("Saved pseudobulk_volcano_plots.pdf")

# =============================================================================
# 5. GOI SUMMARY ACROSS CELL TYPES
# =============================================================================

message("\nBuilding GOI summary...")

goi_summary <- bind_rows(lapply(names(all_de), function(ct) {
  res <- all_de[[ct]]
  res |>
    filter(gene %in% GOI, !is.na(padj)) |>
    mutate(cell_type = ct) |>
    select(gene, cell_type, log2FoldChange, padj)
}))

if (nrow(goi_summary) > 0) {
  goi_summary$sig_label <- ifelse(goi_summary$padj < 0.05, "*", "")
  goi_summary$fc_capped <- pmax(pmin(goi_summary$log2FoldChange, 4), -4)

  p_goi <- ggplot(goi_summary,
                   aes(x = cell_type, y = gene, fill = fc_capped)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_label), size = 5, vjust = 0.75,
              color = "grey20") +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#B71C1C",
                         midpoint = 0,
                         name = "log2FC\n(capped +/-4)") +
    labs(title    = "GOI: WH vs Intact — pseudobulk DESeq2",
         subtitle = "* = padj < 0.05  |  rows = GOI genes, columns = cell types",
         x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7.5),
          axis.text.y  = element_text(size = 8),
          plot.title   = element_text(face = "bold"))

  pdf(file.path(pb_dir, "GOI_pseudobulk_heatmap.pdf"),
      width = max(8, length(unique(goi_summary$cell_type)) * 0.6 + 3),
      height = 7)
  print(p_goi)
  dev.off()
  message("Saved GOI_pseudobulk_heatmap.pdf")
}

# =============================================================================
# SUMMARY
# =============================================================================

message("\n============================================================")
message("PSEUDOBULK DE SUMMARY — WH vs Intact")
message("------------------------------------------------------------")
for (ct in names(all_de)) {
  res <- all_de[[ct]]
  n_sig <- sum(res$padj < 0.05 & !is.na(res$padj))
  n_up  <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  n_dn  <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  message(sprintf("%-35s  DEGs: %4d  (Up: %3d  Down: %3d)", ct, n_sig, n_up, n_dn))
}
message("------------------------------------------------------------")
message("Outputs saved to: ", pb_dir)
message("============================================================")
