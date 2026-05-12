# =============================================================================
# plot_atlas_figures.R
#
# Generates UMAP figures for the cross-species atlas from the saved RDS.
# Run from comparison/ directory: source("R/plot_atlas_figures.R")
# =============================================================================

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

ATLAS_RDS <- "results/cross_species_atlas.rds"
FIG_DIR   <- "figures/atlas"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

message("Loading atlas...")
atlas <- readRDS(ATLAS_RDS)
message(sprintf("  %d cells | %d genes | %d clusters",
                ncol(atlas), nrow(atlas), length(levels(atlas$seurat_clusters))))

# Cluster species composition summary
tab <- table(atlas$seurat_clusters, atlas$species)
pct <- round(100 * prop.table(tab, 1), 1)
cat("\nCluster species composition (%):\n")
print(pct)

# Flag clusters that are >90% one species
dominant <- apply(pct, 1, max)
species_specific <- names(dominant[dominant > 90])
cat("\nClusters >90% one species:", paste(species_specific, collapse = ", "), "\n\n")

# Species colors — consistent with regeneration spectrum plots
species_pal <- c(Axolotl = "#1565C0", Xenopus = "#00897B", Mouse = "#43A047")

# Condition colors
cond_pal <- c(Regen = "#1565C0", NonRegen = "#B71C1C")

pt_size <- pmin(0.5, 30000 / ncol(atlas))

# =============================================================================
# UMAP by species
# =============================================================================
p_species <- DimPlot(atlas, reduction = "umap.harmony",
                     group.by = "species", cols = species_pal,
                     pt.size = pt_size, shuffle = TRUE) +
  ggtitle("Cross-species atlas — by species",
          subtitle = sprintf("%d cells  |  %d genes  |  Harmony (species)",
                             ncol(atlas), nrow(atlas))) +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave(file.path(FIG_DIR, "atlas_umap_species.png"),
       p_species, width = 8, height = 6, dpi = 200)
message("Saved atlas_umap_species.png")

# =============================================================================
# UMAP by cluster (labelled)
# =============================================================================
p_cluster <- DimPlot(atlas, reduction = "umap.harmony",
                     label = TRUE, repel = TRUE, pt.size = pt_size) +
  ggtitle("Cross-species atlas — clusters",
          subtitle = sprintf("%d clusters  |  res 0.5", length(levels(atlas$seurat_clusters))))

ggsave(file.path(FIG_DIR, "atlas_umap_clusters.png"),
       p_cluster, width = 9, height = 6, dpi = 200)
message("Saved atlas_umap_clusters.png")

# =============================================================================
# Species facet (one panel per species, others in grey)
# =============================================================================
p_split <- DimPlot(atlas, reduction = "umap.harmony",
                   group.by = "species", split.by = "species",
                   cols = species_pal, pt.size = pt_size * 0.8,
                   shuffle = TRUE) +
  ggtitle("Cross-species atlas — species split")

ggsave(file.path(FIG_DIR, "atlas_umap_species_split.png"),
       p_split, width = 15, height = 5, dpi = 200)
message("Saved atlas_umap_species_split.png")

# =============================================================================
# Axolotl cell type overlay
# =============================================================================
has_ct <- !is.na(atlas$cell_type) & nzchar(as.character(atlas$cell_type))
if (sum(has_ct) > 500) {
  p_ct <- DimPlot(atlas, reduction = "umap.harmony",
                  cells.highlight = WhichCells(atlas, expression = !is.na(cell_type) & cell_type != ""),
                  group.by = "cell_type",
                  cells = WhichCells(atlas, expression = !is.na(cell_type) & cell_type != ""),
                  pt.size = pt_size * 1.2) +
    ggtitle("Cross-species atlas — axolotl cell type annotations")
  ggsave(file.path(FIG_DIR, "atlas_umap_celltype_axolotl.png"),
         p_ct, width = 10, height = 7, dpi = 200)
  message("Saved atlas_umap_celltype_axolotl.png")
}

# =============================================================================
# By condition (Regen vs NonRegen, where set)
# =============================================================================
cond_cells <- WhichCells(atlas, expression = !is.na(condition))
if (length(cond_cells) > 100) {
  avail_conds <- intersect(names(cond_pal), unique(atlas$condition[cond_cells]))
  p_cond <- DimPlot(atlas, reduction = "umap.harmony",
                    group.by = "condition",
                    cells = cond_cells,
                    cols = cond_pal[avail_conds],
                    pt.size = pt_size, shuffle = TRUE) +
    ggtitle("Cross-species atlas — regeneration condition")
  ggsave(file.path(FIG_DIR, "atlas_umap_condition.png"),
         p_cond, width = 8, height = 6, dpi = 200)
  message("Saved atlas_umap_condition.png")
}

# =============================================================================
# Species-specific clusters — what markers define them?
# =============================================================================
message("\nFinding markers for species-specific clusters...")
Idents(atlas) <- atlas$seurat_clusters

# Use only shared cells for FindMarkers to avoid species confounding
markers_list <- lapply(species_specific, function(cl) {
  tryCatch({
    m <- FindMarkers(atlas, ident.1 = cl, min.pct = 0.25,
                     logfc.threshold = 0.5, verbose = FALSE)
    m$cluster <- cl
    m$gene    <- rownames(m)
    head(m[order(-m$avg_log2FC), ], 10)
  }, error = function(e) NULL)
})
names(markers_list) <- species_specific

cat("\nTop markers for species-specific clusters:\n")
for (cl in species_specific) {
  m <- markers_list[[cl]]
  if (!is.null(m)) {
    sp_dominant <- colnames(pct)[which.max(pct[cl, ])]
    cat(sprintf("\nCluster %s (%s-dominant):\n", cl, sp_dominant))
    cat(paste(sprintf("  %s (log2FC=%.2f, pct1=%.2f)",
                      m$gene, m$avg_log2FC, m$pct.1),
              collapse = "\n"), "\n")
  }
}

message("\nDone. Figures in ", FIG_DIR)
