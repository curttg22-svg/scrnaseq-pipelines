# =============================================================================
# 02c — Xenopus individual GOI plots: per-gene UMAP and violin figures
# Input:  results/xen_merged.rds (annotated, from 02_xenopus_annotation.R)
# Output: results/individual_goi/<GENE>_featureplot.pdf
#         results/individual_goi/<GENE>_featureplot_split.pdf
#         results/individual_goi/<GENE>_violin.pdf
#         results/early_timepoints/individual_goi/ (same, early subset only)
#
# For each GOI that is detectable in the object (homeolog-aware lookup):
#   - Full UMAP colored by expression (all timepoints)
#   - Same UMAP split across timepoints side by side
#   - Violin plot by cell type, split by timepoint
#
# Both the full 5-timepoint object and the early (0-14 dpa) subset are
# processed so figures can be directly compared. Run 02_xenopus_annotation.R
# and 02b_xenopus_early_timepoints.R first.
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)

results_dir <- "results"

TIMEPOINT_COLORS <- c(
  "0dpa"     = "#2E7D32",
  "3dpa"     = "#00838F",
  "7-14dpa"  = "#1565C0",
  "14dpa"    = "#6A1B9A",
  "14-52dpa" = "#B71C1C"
)

EARLY_TIMEPOINTS <- c("0dpa", "3dpa", "7-14dpa", "14dpa")

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

# Maps each GOI to the homeolog name(s) found in the object
build_goi_map <- function(symbols, obj) {
  lapply(setNames(symbols, symbols), function(g) {
    find_xen_genes(g, obj)
  })
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

plot_gene <- function(obj, gene_hits, gene_name, out_dir, tp_colors,
                      label = "all timepoints") {

  if (length(gene_hits) == 0) {
    message("  Skipping ", gene_name, " — no homeologs found in object")
    return(invisible(NULL))
  }

  # Use the homeolog with the highest total expression if multiple
  if (length(gene_hits) > 1) {
    mat   <- GetAssayData(obj, layer = "data")
    tots  <- rowSums(mat[gene_hits, , drop = FALSE])
    gene_plot <- gene_hits[which.max(tots)]
    message("  ", gene_name, ": using ", gene_plot,
            " (highest total of ", length(gene_hits), " homeologs)")
  } else {
    gene_plot <- gene_hits
  }

  max_val <- max(GetAssayData(obj, layer = "data")[gene_plot, ])
  col_lim <- c(0, max(max_val, 0.1))

  # --- Feature plot: full UMAP ---
  p_full <- FeaturePlot(
    obj, features = gene_plot, reduction = "umap_harmony",
    pt.size = 0.25, order = TRUE, cols = c("lightgrey", "#B2182B")
  ) +
    ggtitle(paste0(gene_name, " (", gene_plot, ") - ", label)) +
    theme(plot.title = element_text(size = 10))

  pdf(file.path(out_dir, paste0(gene_name, "_featureplot.pdf")),
      width = 6, height = 5)
  print(p_full)
  dev.off()

  # --- Feature plot: split by timepoint ---
  n_tp   <- length(unique(obj$timepoint))
  p_split <- FeaturePlot(
    obj, features = gene_plot, reduction = "umap_harmony",
    split.by = "timepoint", pt.size = 0.2, order = TRUE,
    cols = c("lightgrey", "#B2182B")
  ) &
    theme(plot.title  = element_text(size = 9),
          strip.text  = element_text(size = 9, face = "bold"),
          legend.position = "right")

  pdf(file.path(out_dir, paste0(gene_name, "_featureplot_split.pdf")),
      width = n_tp * 4, height = 4.5)
  print(p_split)
  dev.off()

  # --- Violin: expression by cell type, split by timepoint ---
  p_vln <- VlnPlot(
    obj, features = gene_plot, group.by = "cell_type",
    split.by = "timepoint", split.plot = TRUE,
    pt.size = 0, cols = tp_colors
  ) +
    ggtitle(paste0(gene_name, " - expression by cell type (", label, ")")) +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 7.5),
          legend.text     = element_text(size = 8),
          legend.title    = element_text(size = 8),
          plot.title      = element_text(size = 9))

  pdf(file.path(out_dir, paste0(gene_name, "_violin.pdf")),
      width = 14, height = 5)
  print(p_vln)
  dev.off()

  message("  Saved: ", gene_name)
}

# =============================================================================
# 1. FULL DATASET (all 5 timepoints)
# =============================================================================

message("Loading annotated object...")
xen <- readRDS(file.path(results_dir, "xen_merged.rds"))
Idents(xen) <- xen$cell_type

out_full <- file.path(results_dir, "individual_goi")
dir.create(out_full, showWarnings = FALSE)

goi_map_full <- build_goi_map(GOI, xen)
found_full   <- sum(sapply(goi_map_full, length) > 0)
message("\nFull dataset — GOI with homeologs found: ", found_full, "/", length(GOI))
message("Generating figures for full dataset (", ncol(xen), " cells)...")

tp_colors_full <- TIMEPOINT_COLORS[levels(xen$timepoint)]

for (gene in GOI) {
  plot_gene(xen, goi_map_full[[gene]], gene, out_full,
            tp_colors_full, label = "all timepoints")
}

# =============================================================================
# 2. EARLY SUBSET (0-14 dpa, excluding 14-52dpa)
# =============================================================================

out_early <- file.path(results_dir, "early_timepoints", "individual_goi")
dir.create(out_early, showWarnings = FALSE)

xen_early <- subset(xen, subset = timepoint %in% EARLY_TIMEPOINTS)
xen_early$timepoint <- droplevels(xen_early$timepoint)
xen_early$cell_type <- droplevels(xen_early$cell_type)
Idents(xen_early)   <- xen_early$cell_type

goi_map_early <- build_goi_map(GOI, xen_early)
message("\nEarly subset (", ncol(xen_early), " cells) — generating figures...")

tp_colors_early <- TIMEPOINT_COLORS[levels(xen_early$timepoint)]

for (gene in GOI) {
  plot_gene(xen_early, goi_map_early[[gene]], gene, out_early,
            tp_colors_early, label = "0-14 dpa")
}

message("\nDone.")
message("Full dataset figures -> ", out_full)
message("Early subset figures -> ", out_early)
