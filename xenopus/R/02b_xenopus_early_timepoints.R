# =============================================================================
# 02b — Xenopus early-timepoint subset analysis
# Input:  results/xen_merged.rds (annotated, from 02_xenopus_annotation.R)
# Output: results/early_timepoints/  (canonical dotplot, UMAP, GOI scores)
#
# Subsets the fully annotated object to the four early timepoints:
#   0dpa, 3dpa, 7-14dpa, 14dpa
# Excludes the 14-52dpa late/maturation pool (12,948 cells, ~34% of total).
#
# Cell type labels from 02_xenopus_annotation.R carry over without
# re-clustering — the same 20 cluster assignments apply, and any cell types
# that drop to 0 cells after subsetting are dropped from the factor.
# Run 02_xenopus_annotation.R first if xen_merged.rds does not yet have
# the cell_type metadata column.
# =============================================================================

library(Seurat)
library(ggplot2)

results_dir <- "results"
out_dir     <- file.path(results_dir, "early_timepoints")
dir.create(out_dir, showWarnings = FALSE)

EARLY_TIMEPOINTS <- c("0dpa", "3dpa", "7-14dpa", "14dpa")

# =============================================================================
# 1. LOAD & SUBSET
# =============================================================================

message("Loading annotated object...")
xen <- readRDS(file.path(results_dir, "xen_merged.rds"))

if (!"cell_type" %in% colnames(xen@meta.data))
  stop("cell_type column not found — run 02_xenopus_annotation.R first")

message("Full object: ", ncol(xen), " cells across ",
        length(unique(xen$timepoint)), " timepoints")
print(table(xen$timepoint))

xen_early <- subset(xen, subset = timepoint %in% EARLY_TIMEPOINTS)
xen_early$timepoint <- droplevels(xen_early$timepoint)
xen_early$cell_type <- droplevels(xen_early$cell_type)
Idents(xen_early)   <- xen_early$cell_type

message("\nEarly subset: ", ncol(xen_early), " cells")
print(table(xen_early$timepoint))
message("\nCell type counts (early subset):")
print(sort(table(xen_early$cell_type), decreasing = TRUE))

# =============================================================================
# 2. CANONICAL MARKER DOTPLOT
# =============================================================================

canonical_markers <- c(
  "col1a1.L","col1a1.S","vim.L","vim.S",
  "col2a1.L","col2a1.S","sox9.L","sox9.S",
  "krt8.L","krt8.S","krt18.L","krt18.S","epcam.L","epcam.S",
  "ptprc.L",
  "s100a8.L",
  "cd79a.L",
  "pecam1.L","pecam1.S","cdh5.L",
  "acta2.L","acta2.S","tagln.L","tagln.S",
  "hba1.L","hba1.S",
  "mki67.L","mki67.S","top2a.L",
  "sox10.L","sox10.S","s100b.L"
)
canonical_markers <- intersect(canonical_markers, rownames(xen_early))

pdf(file.path(out_dir, "early_canonical_marker_dotplot.pdf"),
    width = max(12, length(canonical_markers) * 0.5 + 3), height = 9)
print(
  DotPlot(xen_early, features = canonical_markers, group.by = "cell_type",
          dot.scale = 6, col.min = 0) +
    scale_color_gradient(low = "lightgrey", high = "#08306B") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8)) +
    ggtitle("Xenopus early timepoints (0-14 dpa) - canonical markers by cell type")
)
dev.off()
message("Saved early_canonical_marker_dotplot.pdf")

# =============================================================================
# 3. ANNOTATED UMAP — full and split by timepoint
# =============================================================================

pdf(file.path(out_dir, "early_umap_annotated.pdf"), width = 10, height = 7)
print(
  DimPlot(xen_early, reduction = "umap_harmony", group.by = "cell_type",
          label = TRUE, repel = TRUE, pt.size = 0.3) +
    ggtitle("Xenopus early timepoints (0-14 dpa) - annotated") +
    theme(legend.text = element_text(size = 8))
)
dev.off()

pdf(file.path(out_dir, "early_umap_by_timepoint.pdf"), width = 12, height = 4)
print(
  DimPlot(xen_early, reduction = "umap_harmony", group.by = "timepoint",
          pt.size = 0.3,
          cols = c("0dpa"    = "#2E7D32",
                   "3dpa"    = "#00838F",
                   "7-14dpa" = "#1565C0",
                   "14dpa"   = "#6A1B9A")) +
    ggtitle("Xenopus early timepoints - by timepoint")
)
dev.off()

pdf(file.path(out_dir, "early_umap_split_by_timepoint.pdf"),
    width = 20, height = 5)
print(
  DimPlot(xen_early, reduction = "umap_harmony", group.by = "cell_type",
          split.by = "timepoint", label = TRUE, repel = TRUE,
          pt.size = 0.2, label.size = 2.5) +
    theme(legend.position = "none",
          strip.text      = element_text(size = 10, face = "bold"))
)
dev.off()
message("Saved UMAP figures")

# =============================================================================
# 4. GOI MODULE SCORE
# =============================================================================

GOI <- c(
  "GLI1","PTCH1","PTCH2","SHH","IHH","DHH",
  "GDF5","MSX1","SALL1","GREM1",
  "FREM1","CYTL1","CLDN5","IL11","FGFBP1",
  "C3","LDLR","ALOX15B","MSMO1","CYP51A1","FDPS"
)

find_xen_genes <- function(symbols, obj) {
  all_genes <- rownames(obj)
  found <- c()
  for (g in symbols) {
    hits <- grep(paste0("^", g, "(\\.L|\\.S|\\.L\\.[0-9]+|\\.S\\.[0-9]+)?$"),
                 all_genes, ignore.case = TRUE, value = TRUE)
    found <- c(found, hits)
  }
  unique(found)
}

xen_goi <- find_xen_genes(GOI, xen_early)
message("GOI homeologs found: ", length(xen_goi))

if (length(xen_goi) > 0) {
  xen_early <- AddModuleScore(xen_early, features = list(xen_goi),
                               name = "GOI_score")

  pdf(file.path(out_dir, "early_GOI_module_score_UMAP.pdf"), width = 8, height = 6)
  print(
    FeaturePlot(xen_early, features = "GOI_score1",
                reduction = "umap_harmony", pt.size = 0.3, order = TRUE) +
      scale_color_gradient2(low = "lightgrey", mid = "#FDAE61",
                            high = "#B2182B", midpoint = 0) +
      ggtitle("Xenopus early timepoints - GOI module score") +
      theme(legend.title = element_text(size = 8))
  )
  dev.off()

  pdf(file.path(out_dir, "early_GOI_module_score_violin.pdf"),
      width = 12, height = 5)
  print(
    VlnPlot(xen_early, features = "GOI_score1", group.by = "cell_type",
            pt.size = 0) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "none") +
      ggtitle("GOI module score by cell type - early timepoints (0-14 dpa)")
  )
  dev.off()

  # Split UMAP: one panel per timepoint, GOI score as color
  pdf(file.path(out_dir, "early_GOI_score_split_by_timepoint.pdf"),
      width = 20, height = 5)
  print(
    FeaturePlot(xen_early, features = "GOI_score1",
                reduction = "umap_harmony", split.by = "timepoint",
                pt.size = 0.2, order = TRUE) &
      scale_color_gradient2(low = "lightgrey", mid = "#FDAE61",
                            high = "#B2182B", midpoint = 0) &
      theme(legend.position = "right")
  )
  dev.off()

  message("GOI module score figures saved")
}

message("\nAll early-timepoint figures saved to: ", out_dir)
