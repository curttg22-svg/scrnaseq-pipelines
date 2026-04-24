# =============================================================================
# 07 — Exploratory DE: EB vs Intact and MB vs Intact (per cell type)
# Input:  results/axo_integrated.rds
# Output: results/exploratory_de/  — CSV tables + PDF figures
#
# IMPORTANT — STATISTICAL CAVEAT:
#   EB (SCP499) and MB (SCP500) each have only ONE biological sample.
#   Pseudobulk DESeq2 is NOT valid without replicates.
#   This script uses cell-level Wilcoxon rank-sum test per harmonized cell type
#   as an EXPLORATORY analysis only. Results should NOT be interpreted as
#   statistically confirmed differential expression.
#
#   Stricter thresholds are applied to compensate for inflated significance:
#     - |log2FC| > 1 (2-fold change minimum)
#     - p_val_adj < 0.01 (more stringent than standard 0.05)
#     - min.pct = 0.2 (gene detected in at least 20% of cells in one group)
#
#   Use these results to generate hypotheses and prioritize genes for
#   validation — not as standalone evidence of differential expression.
#   Confirmatory analysis would require additional biological replicates.
#
# CIRBP NOTE:
#   As observed in the WH pseudobulk analysis, CIRBP is likely a global
#   wound/stress response gene or batch artifact. Its appearance here
#   should be interpreted with the same caution.
#
# Reference: Squair et al. 2021, Nat Commun — single-sample comparisons
#   are fundamentally limited; treat as exploratory.
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

results_dir <- "results"
exp_dir     <- file.path(results_dir, "exploratory_de")
dir.create(exp_dir, showWarnings = FALSE, recursive = TRUE)

WH_REPS <- c("WH_N4", "WH_N5", "WH_N6")

GOI <- c(
  "GLI1","PTCH1","PTCH2","SHH","IHH","DHH",
  "GDF5","MSX1","SALL1","GREM1",
  "FREM1","CYTL1","CLDN5","IL11","FGFBP1",
  "C3","LDLR","ALOX15B","MSMO1","CYP51A1","FDPS"
)

# Same harmonization map as script 06
CT_HARMONIZE <- data.frame(
  intact_label = c(
    "SSCs","Fibroblast","Tcells","Early_B","RBCs","Myeloid",
    "BasalEpidermis","Intermediate_epidermis","Proliferating_Epidermis",
    "Epidermal_Langerhans"
  ),
  eb_mb_label = c(
    "SSC","Fibroblast-like_blastema_","T_cell","Early_B_cell","Erythrocyte",
    "Recruited_Macrophage","Basal WE","Intermediate_WE_","Intermediate_WE_",
    "Dendritic_Cell"
  ),
  harmonized = c(
    "SSC_fibroblast","Fibroblast_blastema","T_cell","Early_B_cell",
    "Erythrocyte_RBC","Myeloid_macrophage","Basal_epidermis",
    "Intermediate_epidermis","Intermediate_epidermis","Dendritic_Langerhans"
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# 1. LOAD & PREPARE
# =============================================================================

message("Loading axo_integrated.rds...")
axo <- readRDS(file.path(results_dir, "axo_integrated.rds"))
axo <- JoinLayers(axo)

# Assign harmonized cell type labels
axo$ct_harmonized <- NA_character_

# Intact labels
for (i in seq_len(nrow(CT_HARMONIZE))) {
  idx <- axo$timepoint == "Intact" &
         axo$paper_cluster == CT_HARMONIZE$intact_label[i]
  axo$ct_harmonized[idx] <- CT_HARMONIZE$harmonized[i]
}
# EB labels
for (i in seq_len(nrow(CT_HARMONIZE))) {
  idx <- axo$timepoint == "EB" &
         axo$paper_cluster == CT_HARMONIZE$eb_mb_label[i]
  axo$ct_harmonized[idx] <- CT_HARMONIZE$harmonized[i]
}
# MB labels
for (i in seq_len(nrow(CT_HARMONIZE))) {
  idx <- axo$timepoint == "MB" &
         axo$paper_cluster == CT_HARMONIZE$eb_mb_label[i]
  axo$ct_harmonized[idx] <- CT_HARMONIZE$harmonized[i]
}

message("Harmonized cell type counts per timepoint:")
print(table(axo$ct_harmonized, axo$timepoint)[,
      c("Intact","EB","MB"), drop = FALSE])

# =============================================================================
# 2. RUN WILCOXON PER CELL TYPE — EB vs INTACT, MB vs INTACT
# =============================================================================

run_exploratory_de <- function(axo, comparison, tp, min_cells = 20) {
  message("\n--- ", comparison, " ---")

  # Subset to Intact + target timepoint, harmonized cells only
  axo_sub <- subset(axo,
    timepoint %in% c("Intact", tp) & !is.na(ct_harmonized))

  Idents(axo_sub) <- axo_sub$ct_harmonized
  ct_levels <- unique(axo_sub$ct_harmonized)

  all_de <- list()

  for (ct in ct_levels) {
    ct_obj <- subset(axo_sub, idents = ct)
    Idents(ct_obj) <- ct_obj$timepoint

    n_intact <- sum(ct_obj$timepoint == "Intact")
    n_tp     <- sum(ct_obj$timepoint == tp)

    if (n_intact < min_cells || n_tp < min_cells) {
      message("  Skipping ", ct, " (Intact: ", n_intact, ", ", tp, ": ", n_tp, ")")
      next
    }

    message("  ", ct, " (Intact: ", n_intact, ", ", tp, ": ", n_tp, ")")

    de <- tryCatch(
      FindMarkers(
        ct_obj,
        ident.1         = tp,
        ident.2         = "Intact",
        test.use        = "wilcox",
        min.pct         = 0.2,
        logfc.threshold = 1.0,   # stricter: 2-fold minimum
        verbose         = FALSE
      ),
      error = function(e) { message("    Error: ", e$message); NULL }
    )
    if (is.null(de)) next

    de$gene       <- rownames(de)
    de$cell_type  <- ct
    de$comparison <- comparison

    # Apply stricter p-value threshold
    n_sig <- sum(de$p_val_adj < 0.01, na.rm = TRUE)
    message("    DEGs (p_adj<0.01, |log2FC|>1): ", n_sig)

    all_de[[ct]] <- de
    ct_label <- gsub("[^A-Za-z0-9_]", "_", ct)
    write.csv(de,
              file.path(exp_dir,
                        paste0(ct_label, "_", comparison, "_exploratory.csv")),
              row.names = FALSE)
  }

  all_de
}

eb_de <- run_exploratory_de(axo, "EB_vs_Intact", "EB")
mb_de <- run_exploratory_de(axo, "MB_vs_Intact", "MB")

# =============================================================================
# 3. VOLCANO PLOTS
# =============================================================================

make_exp_volcano <- function(de, title, goi = GOI) {
  de <- de[!is.na(de$p_val_adj), ]
  de$sig <- "NS"
  de$sig[de$p_val_adj < 0.01 & de$avg_log2FC >  1] <- "Up"
  de$sig[de$p_val_adj < 0.01 & de$avg_log2FC < -1] <- "Down"
  de$sig <- factor(de$sig, levels = c("Up","Down","NS"))

  de$log10p  <- pmin(-log10(de$p_val_adj + 1e-300), 30)
  de$fc_plot <- pmax(pmin(de$avg_log2FC, 8), -8)

  de_named  <- de |> filter(!grepl("^c[0-9]|-g[0-9]|-i[0-9]|\\^", gene))
  top_up    <- de_named |> filter(sig == "Up")   |> arrange(desc(avg_log2FC)) |> head(8)
  top_down  <- de_named |> filter(sig == "Down")  |> arrange(avg_log2FC)      |> head(8)
  goi_hits  <- de       |> filter(gene %in% goi, p_val_adj < 0.01)
  label_df  <- bind_rows(top_up, top_down, goi_hits) |> distinct(gene, .keep_all = TRUE)

  ggplot(de, aes(x = fc_plot, y = log10p, color = sig)) +
    geom_point(size = 0.9, alpha = 0.5) +
    geom_point(data = label_df, size = 2, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = label_df, aes(label = gene),
      size = 2.8, max.overlaps = 15,
      box.padding = 0.4, color = "grey20"
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    scale_color_manual(
      values = c("Up" = "#B71C1C", "Down" = "#1565C0", "NS" = "grey75"),
      name = NULL
    ) +
    labs(title    = title,
         subtitle = paste0("EXPLORATORY — single sample, no replicates",
                           "\nUp: ", sum(de$sig == "Up"),
                           "  |  Down: ", sum(de$sig == "Down"),
                           "  |  p_adj<0.01, |log2FC|>1"),
         x = "log2 fold-change vs Intact (capped at +/-8)",
         y = "-log10 adjusted p-value (capped at 30)") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(color = "#B71C1C", size = 8),
          legend.position = "right")
}

message("\nGenerating volcano plots...")

pdf(file.path(exp_dir, "EB_vs_Intact_volcano_plots.pdf"), width = 7, height = 5.5)
for (ct in names(eb_de)) {
  print(make_exp_volcano(eb_de[[ct]],
                         title = paste0(ct, " - EB vs Intact (exploratory)")))
}
dev.off()
message("Saved EB_vs_Intact_volcano_plots.pdf")

pdf(file.path(exp_dir, "MB_vs_Intact_volcano_plots.pdf"), width = 7, height = 5.5)
for (ct in names(mb_de)) {
  print(make_exp_volcano(mb_de[[ct]],
                         title = paste0(ct, " - MB vs Intact (exploratory)")))
}
dev.off()
message("Saved MB_vs_Intact_volcano_plots.pdf")

# =============================================================================
# 4. GOI HEATMAPS — EB and MB side by side
# =============================================================================

message("\nBuilding GOI heatmaps...")

build_goi_df <- function(de_list, comparison) {
  bind_rows(lapply(names(de_list), function(ct) {
    de_list[[ct]] |>
      filter(gene %in% GOI, !is.na(p_val_adj)) |>
      mutate(cell_type  = ct,
             comparison = comparison,
             sig_label  = ifelse(p_val_adj < 0.01, "*", ""),
             fc_capped  = pmax(pmin(avg_log2FC, 4), -4)) |>
      select(gene, cell_type, comparison, avg_log2FC, p_val_adj,
             sig_label, fc_capped)
  }))
}

eb_goi <- build_goi_df(eb_de, "EB vs Intact")
mb_goi <- build_goi_df(mb_de, "MB vs Intact")
goi_combined <- bind_rows(eb_goi, mb_goi)

make_goi_heatmap <- function(df, title) {
  if (nrow(df) == 0) {
    message("  No GOI hits for: ", title)
    return(NULL)
  }
  ggplot(df, aes(x = cell_type, y = gene, fill = fc_capped)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_label), size = 5, vjust = 0.75,
              color = "grey20") +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#B71C1C",
                         midpoint = 0, name = "log2FC\n(capped +/-4)") +
    labs(title    = title,
         subtitle = "* = p_adj < 0.01  |  EXPLORATORY: single sample, no replicates",
         x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7.5),
          axis.text.y  = element_text(size = 8),
          plot.title   = element_text(face = "bold"),
          plot.subtitle = element_text(color = "#B71C1C", size = 7.5))
}

p_eb <- make_goi_heatmap(eb_goi, "GOI: EB vs Intact (exploratory)")
p_mb <- make_goi_heatmap(mb_goi, "GOI: MB vs Intact (exploratory)")

if (!is.null(p_eb)) {
  pdf(file.path(exp_dir, "GOI_EB_vs_Intact_heatmap.pdf"),
      width = max(8, length(unique(eb_goi$cell_type)) * 0.7 + 3), height = 7)
  print(p_eb)
  dev.off()
  message("Saved GOI_EB_vs_Intact_heatmap.pdf")
}

if (!is.null(p_mb)) {
  pdf(file.path(exp_dir, "GOI_MB_vs_Intact_heatmap.pdf"),
      width = max(8, length(unique(mb_goi$cell_type)) * 0.7 + 3), height = 7)
  print(p_mb)
  dev.off()
  message("Saved GOI_MB_vs_Intact_heatmap.pdf")
}

# =============================================================================
# SUMMARY
# =============================================================================

message("\n============================================================")
message("EXPLORATORY DE SUMMARY (single-sample, treat as hypothesis-generating)")
message("------------------------------------------------------------")
message("EB vs Intact:")
for (ct in names(eb_de)) {
  n <- sum(eb_de[[ct]]$p_val_adj < 0.01, na.rm = TRUE)
  message(sprintf("  %-30s  DEGs: %d", ct, n))
}
message("MB vs Intact:")
for (ct in names(mb_de)) {
  n <- sum(mb_de[[ct]]$p_val_adj < 0.01, na.rm = TRUE)
  message(sprintf("  %-30s  DEGs: %d", ct, n))
}
message("------------------------------------------------------------")
message("Outputs saved to: ", exp_dir)
message("============================================================")
