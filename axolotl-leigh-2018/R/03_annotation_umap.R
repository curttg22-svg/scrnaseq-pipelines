# =============================================================================
# 03 — UMAPs, composition, GOI dotplots / violins / featureplots
# Input:  results/axo_integrated.rds
# Output: PDF figures in results/
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

results_dir <- "results"

TIMEPOINT_ORDER <- c("Intact", "WH_N4", "WH_N5", "WH_N6", "EB", "MB")
TIMEPOINT_COLORS <- c(
  "Intact" = "#2E7D32",
  "WH_N4"  = "#00838F",
  "WH_N5"  = "#006064",
  "WH_N6"  = "#004D40",
  "EB"     = "#1565C0",
  "MB"     = "#6A1B9A"
)

GOI <- c(
  "GLI1","PTCH1","PTCH2","SHH","IHH","DHH",
  "GDF5","MSX1","SALL1","GREM1",
  "FREM1","CYTL1","CLDN5","IL11","FGFBP1",
  "C3","LDLR","ALOX15B","MSMO1","CYP51A1","FDPS"
)

axo <- readRDS(file.path(results_dir, "axo_integrated.rds"))
DefaultDimReduc(axo) <- "umap_harmony"

# =============================================================================
# UMAPs
# =============================================================================

p_time <- DimPlot(axo, group.by = "timepoint", pt.size = 0.3,
                   cols = TIMEPOINT_COLORS, shuffle = TRUE) +
  ggtitle("Axolotl regeneration time course\nHarmony-integrated UMAP") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))

p_ct <- DimPlot(axo, group.by = "paper_cluster", pt.size = 0.3,
                 label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("Paper cell type annotations") +
  theme_classic(base_size = 11) +
  theme(plot.title   = element_text(hjust = 0.5),
        legend.text  = element_text(size = 7))

p_study <- DimPlot(axo, group.by = "study_id", pt.size = 0.3, shuffle = TRUE) +
  ggtitle("Study of origin (batch check)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file.path(results_dir, "axo_UMAP_overview.pdf"), width = 16, height = 6)
print(p_time + p_ct + p_study + plot_layout(ncol = 3))
dev.off()
message("Saved axo_UMAP_overview.pdf")

# Per-timepoint paper UMAPs
make_paper_umap <- function(obj, tp, color) {
  sub_obj <- subset(obj, timepoint == tp)
  DimPlot(sub_obj, reduction = "paper_umap", group.by = "paper_cluster",
          pt.size = 0.4, label = TRUE, label.size = 2.5, repel = TRUE) +
    ggtitle(tp) +
    theme_classic(base_size = 9) +
    theme(plot.title      = element_text(hjust = 0.5, color = color),
          legend.text     = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"))
}

pdf(file.path(results_dir, "axo_paper_UMAPs.pdf"), width = 18, height = 10)
p_intact <- make_paper_umap(axo, "Intact", TIMEPOINT_COLORS["Intact"])
wh_sub   <- subset(axo, study_id == "SCP489")
p_wh <- DimPlot(wh_sub, reduction = "paper_umap", group.by = "paper_cluster",
                 pt.size = 0.4, label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("Wound Healing (N4 + N5 + N6)") +
  theme_classic(base_size = 9) +
  theme(plot.title  = element_text(hjust = 0.5),
        legend.text = element_text(size = 6))
p_eb <- make_paper_umap(axo, "EB", TIMEPOINT_COLORS["EB"])
p_mb <- make_paper_umap(axo, "MB", TIMEPOINT_COLORS["MB"])
print((p_intact | p_wh) / (p_eb | p_mb))
dev.off()
message("Saved axo_paper_UMAPs.pdf")

pdf(file.path(results_dir, "axo_UMAP_by_timepoint.pdf"), width = 20, height = 4)
print(
  DimPlot(axo, group.by = "paper_cluster", split.by = "timepoint",
          ncol = 6, pt.size = 0.2, label = FALSE) +
    theme_classic(base_size = 8) +
    theme(legend.position = "none")
)
dev.off()
message("Saved axo_UMAP_by_timepoint.pdf")

# =============================================================================
# Composition bar chart
# =============================================================================

comp_df <- as.data.frame(table(axo$timepoint, axo$paper_cluster))
colnames(comp_df) <- c("Timepoint", "CellType", "Count")
comp_df <- comp_df |>
  group_by(Timepoint) |>
  mutate(Proportion = Count / sum(Count)) |>
  ungroup()

p_comp <- ggplot(comp_df, aes(x = Timepoint, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(limits = TIMEPOINT_ORDER) +
  theme_classic(base_size = 11) +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        legend.text      = element_text(size = 7),
        legend.key.size  = unit(0.3, "cm"),
        plot.title       = element_text(hjust = 0.5)) +
  ggtitle("Cell type composition across timepoints") +
  ylab("Proportion") + xlab("")

pdf(file.path(results_dir, "axo_composition_barplot.pdf"), width = 10, height = 6)
print(p_comp)
dev.off()
message("Saved axo_composition_barplot.pdf")

# =============================================================================
# GOI dotplots, featureplots, violins
# =============================================================================

goi_present <- intersect(GOI, rownames(axo))
message("GOI in object: ", length(goi_present), "/", length(GOI))

if (length(goi_present) > 0) {

  pdf(file.path(results_dir, "axo_GOI_dotplot_celltypes.pdf"),
      width = max(10, length(goi_present) * 0.7 + 3), height = 12)
  print(
    DotPlot(axo, features = goi_present, group.by = "paper_cluster",
            dot.scale = 6, col.min = 0) +
      scale_color_gradient(low = "lightgrey", high = "#08306B") +
      ggtitle("GOI expression across cell types (all timepoints)") +
      theme_classic(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_text(size = 8),
            plot.title  = element_text(hjust = 0.5))
  )
  dev.off()
  message("Saved axo_GOI_dotplot_celltypes.pdf")

  pdf(file.path(results_dir, "axo_GOI_dotplot_timepoints.pdf"),
      width = max(10, length(goi_present) * 0.7 + 3), height = 5)
  print(
    DotPlot(axo, features = goi_present, group.by = "timepoint",
            dot.scale = 7, col.min = 0) +
      scale_color_gradient(low = "lightgrey", high = "#B2182B") +
      ggtitle("GOI expression across timepoints") +
      theme_classic(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            plot.title  = element_text(hjust = 0.5))
  )
  dev.off()
  message("Saved axo_GOI_dotplot_timepoints.pdf")

  batches <- split(goi_present, ceiling(seq_along(goi_present) / 6))
  pdf(file.path(results_dir, "axo_GOI_featureplots.pdf"), width = 14, height = 10)
  for (b in batches) {
    print(
      FeaturePlot(axo, features = b, ncol = 3, pt.size = 0.2,
                  order = TRUE, cols = c("lightgrey", "darkblue")) &
        theme_classic(base_size = 9)
    )
  }
  dev.off()
  message("Saved axo_GOI_featureplots.pdf")

  pdf(file.path(results_dir, "axo_GOI_violin_timepoints.pdf"),
      width = 14, height = ceiling(length(goi_present) / 3) * 3.5)
  print(
    VlnPlot(axo, features = goi_present, group.by = "timepoint",
            pt.size = 0, cols = TIMEPOINT_COLORS, ncol = 3) &
      theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 8),
            axis.title     = element_blank(),
            legend.position = "none")
  )
  dev.off()
  message("Saved axo_GOI_violin_timepoints.pdf")
}

message("\nAll figures saved to ", results_dir)
