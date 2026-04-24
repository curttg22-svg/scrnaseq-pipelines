# =============================================================================
# 05 -Differential expression across regeneration timepoints
# Input:  results/axo_integrated.rds
# Output: results/de/  -CSV tables + PDF figures per comparison
#
# Compares WH, EB, and MB each against Intact as baseline.
# Also runs within-cell-type DE for the most abundant cell populations.
#
# Adapted from NXR-2026 Exercise 07 (marker method comparison) applied to
# the Leigh 2018 multi-timepoint axolotl regeneration dataset.
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

results_dir <- "results"
de_dir      <- file.path(results_dir, "de")
dir.create(de_dir, showWarnings = FALSE, recursive = TRUE)

TIMEPOINT_ORDER <- c("Intact", "WH", "EB", "MB")
WH_REPS         <- c("WH_N4", "WH_N5", "WH_N6")

TIMEPOINT_COLORS <- c(
  "Intact" = "#2E7D32",
  "WH"     = "#00838F",
  "EB"     = "#1565C0",
  "MB"     = "#6A1B9A"
)

GOI <- c(
  "GLI1","PTCH1","PTCH2","SHH","IHH","DHH",
  "GDF5","MSX1","SALL1","GREM1",
  "FREM1","CYTL1","CLDN5","IL11","FGFBP1",
  "C3","LDLR","ALOX15B","MSMO1","CYP51A1","FDPS"
)

# =============================================================================
# 1. LOAD & PREPARE
# =============================================================================

message("Loading axo_integrated.rds...")
axo <- readRDS(file.path(results_dir, "axo_integrated.rds"))
axo <- JoinLayers(axo)

# Pool WH replicates into single "WH" stage
axo$timepoint_pooled <- as.character(axo$timepoint)
axo$timepoint_pooled[axo$timepoint_pooled %in% WH_REPS] <- "WH"
axo$timepoint_pooled <- factor(axo$timepoint_pooled, levels = TIMEPOINT_ORDER)

message("Cells per timepoint:")
print(table(axo$timepoint_pooled))

# =============================================================================
# 2. GLOBAL DE: each timepoint vs Intact
# =============================================================================

message("\nRunning global DE (each timepoint vs Intact)...")
Idents(axo) <- axo$timepoint_pooled

comparisons <- c("WH", "EB", "MB")
de_results  <- list()

for (tp in comparisons) {
  message("  ", tp, " vs Intact...")
  de <- FindMarkers(
    axo,
    ident.1         = tp,
    ident.2         = "Intact",
    test.use        = "wilcox",
    min.pct         = 0.1,
    logfc.threshold = 0.25,
    verbose         = FALSE
  )
  de$gene       <- rownames(de)
  de$comparison <- paste0(tp, "_vs_Intact")
  de_results[[tp]] <- de

  write.csv(de,
            file.path(de_dir, paste0(tp, "_vs_Intact_DE.csv")),
            row.names = FALSE)
  message("    ", nrow(de), " DEGs (p_val_adj < 0.05: ",
          sum(de$p_val_adj < 0.05, na.rm = TRUE), ")")
}

# =============================================================================
# 3. VOLCANO PLOTS
# =============================================================================

message("\nGenerating volcano plots...")

make_volcano <- function(de, title, goi = GOI) {
  de$sig <- "NS"
  de$sig[de$p_val_adj < 0.05 & de$avg_log2FC >  0.5] <- "Up"
  de$sig[de$p_val_adj < 0.05 & de$avg_log2FC < -0.5] <- "Down"
  de$sig <- factor(de$sig, levels = c("Up","Down","NS"))

  # Only label clean gene symbols — exclude Trinity IDs (contain digits after dash)
  de_named <- de |> filter(!grepl("^c[0-9]|-g[0-9]|-i[0-9]|\\^", gene))

  # Top 8 up + top 8 down by fold-change among named genes + significant GOI
  top_up   <- de_named |> filter(sig == "Up")   |> arrange(desc(avg_log2FC)) |> head(8)
  top_down <- de_named |> filter(sig == "Down")  |> arrange(avg_log2FC)      |> head(8)
  goi_hits <- de       |> filter(gene %in% goi, p_val_adj < 0.05)
  label_df <- bind_rows(top_up, top_down, goi_hits) |> distinct(gene, .keep_all = TRUE)

  ggplot(de, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300),
                  color = sig)) +
    geom_point(size = 0.8, alpha = 0.5) +
    geom_point(data = label_df, size = 1.5, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = label_df, aes(label = gene),
      size = 2.8, max.overlaps = 20,
      box.padding = 0.3, color = "grey20"
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    scale_color_manual(values = c("Up" = "#B71C1C", "Down" = "#1565C0",
                                   "NS" = "grey75"),
                       name = NULL) +
    labs(title    = title,
         subtitle = paste0("Up: ", sum(de$sig == "Up"),
                           "  |  Down: ", sum(de$sig == "Down"),
                           "  |  NS: ", sum(de$sig == "NS")),
         x = "Average log2 fold-change",
         y = "-log10 adjusted p-value") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          legend.position = "right")
}

# Check ggrepel is available
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  message("  ggrepel not installed -installing...")
  install.packages("ggrepel")
}
library(ggrepel)

pdf(file.path(de_dir, "volcano_plots.pdf"), width = 7, height = 5)
for (tp in comparisons) {
  print(make_volcano(de_results[[tp]],
                     title = paste0(tp, " vs Intact -all cell types")))
}
dev.off()
message("Saved de/volcano_plots.pdf")

# =============================================================================
# 4. TOP DEG HEATMAP
# =============================================================================

message("\nGenerating top DEG heatmap...")

# Top 15 up + 15 down per comparison
top_genes <- bind_rows(lapply(comparisons, function(tp) {
  de <- de_results[[tp]] |> filter(p_val_adj < 0.05)
  up   <- de |> arrange(desc(avg_log2FC)) |> head(15)
  down <- de |> arrange(avg_log2FC)        |> head(15)
  bind_rows(up, down)
})) |>
  distinct(gene) |>
  pull(gene)

top_genes <- intersect(top_genes, rownames(axo))
message("  Unique top genes for heatmap: ", length(top_genes))

if (length(top_genes) > 0) {
  # Scale just the genes needed -DoHeatmap requires scale.data
  axo <- ScaleData(axo, features = top_genes, verbose = FALSE)

  pdf(file.path(de_dir, "top_deg_heatmap.pdf"), width = 10, height = 8)
  print(
    DoHeatmap(axo,
              features = top_genes,
              group.by = "timepoint_pooled",
              group.colors = TIMEPOINT_COLORS,
              size = 3, angle = 0) +
      scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#B71C1C",
                            midpoint = 0, name = "Scaled\nexpr.") +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
  message("Saved de/top_deg_heatmap.pdf")
}

# =============================================================================
# 5. GOI EXPRESSION CHANGE ACROSS TIMEPOINTS
# =============================================================================

message("\nGOI fold-change summary...")

goi_present <- intersect(GOI, rownames(axo))

goi_fc <- bind_rows(lapply(comparisons, function(tp) {
  de <- de_results[[tp]]
  de |>
    filter(gene %in% goi_present) |>
    mutate(timepoint = tp) |>
    select(gene, timepoint, avg_log2FC, p_val_adj)
}))

goi_fc$timepoint <- factor(goi_fc$timepoint, levels = comparisons)
goi_fc$sig_label <- ifelse(goi_fc$p_val_adj < 0.05, "*", "")

p_goi <- ggplot(goi_fc, aes(x = timepoint, y = gene, fill = avg_log2FC)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 5, vjust = 0.75, color = "grey20") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#B71C1C",
                       midpoint = 0, name = "log2FC\nvs Intact") +
  labs(title    = "GOI expression change vs Intact",
       subtitle = "* = p_adj < 0.05 | X-axis: timepoint vs Intact",
       x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title   = element_text(face = "bold"),
        axis.text.x  = element_text(size = 9),
        axis.text.y  = element_text(size = 8))

pdf(file.path(de_dir, "GOI_foldchange_heatmap.pdf"), width = 6, height = 7)
print(p_goi)
dev.off()
message("Saved de/GOI_foldchange_heatmap.pdf")

# =============================================================================
# 6. WITHIN-CELL-TYPE DE (top 3 most abundant cell types)
# =============================================================================

message("\nRunning within-cell-type DE...")

ct_counts  <- sort(table(axo$paper_cluster), decreasing = TRUE)
top_3_ct   <- names(ct_counts)[1:3]
message("  Top 3 cell types: ", paste(top_3_ct, collapse = ", "))

Idents(axo) <- axo$paper_cluster

for (ct in top_3_ct) {
  ct_obj <- subset(axo, idents = ct)
  Idents(ct_obj) <- ct_obj$timepoint_pooled

  ct_label <- gsub("[^A-Za-z0-9_]", "_", ct)
  ct_de_list <- list()

  for (tp in comparisons) {
    n_tp     <- sum(ct_obj$timepoint_pooled == tp)
    n_intact <- sum(ct_obj$timepoint_pooled == "Intact")
    if (n_tp < 10 || n_intact < 10) {
      message("  Skipping ", ct, " -", tp, " (too few cells: ",
              n_tp, " vs ", n_intact, ")")
      next
    }
    de <- tryCatch(
      FindMarkers(ct_obj, ident.1 = tp, ident.2 = "Intact",
                  test.use = "wilcox", min.pct = 0.1,
                  logfc.threshold = 0.25, verbose = FALSE),
      error = function(e) { message("  Error: ", e$message); NULL }
    )
    if (is.null(de)) next
    de$gene       <- rownames(de)
    de$comparison <- paste0(tp, "_vs_Intact")
    ct_de_list[[tp]] <- de
    write.csv(de,
              file.path(de_dir, paste0(ct_label, "_", tp, "_vs_Intact_DE.csv")),
              row.names = FALSE)
    message("  ", ct, " | ", tp, " vs Intact: ", nrow(de), " DEGs")
  }

  if (length(ct_de_list) == 0) next

  # Volcano for this cell type
  pdf(file.path(de_dir, paste0(ct_label, "_volcano_plots.pdf")),
      width = 7, height = 5)
  for (tp in names(ct_de_list)) {
    print(make_volcano(ct_de_list[[tp]],
                       title = paste0(ct, " -", tp, " vs Intact")))
  }
  dev.off()
  message("  Saved: ", ct_label, "_volcano_plots.pdf")
}

# =============================================================================
# SUMMARY
# =============================================================================

message("\n============================================================")
message("TIMEPOINT DE SUMMARY")
message("------------------------------------------------------------")
for (tp in comparisons) {
  de <- de_results[[tp]]
  message(tp, " vs Intact | Total DEGs: ", nrow(de),
          " | Sig (p_adj<0.05): ", sum(de$p_val_adj < 0.05, na.rm = TRUE),
          " | Up: ", sum(de$p_val_adj < 0.05 & de$avg_log2FC > 0, na.rm = TRUE),
          " | Down: ", sum(de$p_val_adj < 0.05 & de$avg_log2FC < 0, na.rm = TRUE))
}
message("------------------------------------------------------------")
message("GOI detected in DE: ", length(goi_present), "/", length(GOI))
message("Outputs saved to: ", de_dir)
message("============================================================")
