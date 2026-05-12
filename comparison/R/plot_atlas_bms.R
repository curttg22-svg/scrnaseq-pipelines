# =============================================================================
# plot_atlas_bms.R
#
# 1. Annotates all atlas clusters with cell type labels (species-specific
#    clusters from marker analysis; mixed clusters from majority Xenopus/
#    axolotl label).
# 2. Saves annotated cell_type_atlas column back to the atlas RDS.
# 3. Projects BMS bulk RNA-seq onto the atlas via per-cluster pseudobulk
#    Spearman correlation and generates UMAP figures.
#
# Run from comparison/ directory: source("R/plot_atlas_bms.R")
# =============================================================================

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

ATLAS_RDS  <- "results/cross_species_atlas.rds"
CACHE_RDS  <- "results/pseudobulk_cache.rds"
FIG_DIR    <- "figures/atlas"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load atlas
# =============================================================================

message("Loading atlas...")
atlas <- readRDS(ATLAS_RDS)
message(sprintf("  %d cells | %d genes | %d clusters",
                ncol(atlas), nrow(atlas), length(levels(atlas$seurat_clusters))))

# =============================================================================
# 2. Assign cell type labels
# =============================================================================

# Manual labels for species-specific clusters (from marker analysis)
species_specific_labels <- c(
  "0"  = "Injury CT fibroblast [Mouse]",
  "2"  = "Erythroid cells [Axolotl]",
  "3"  = "Resting digit CT fibroblast [Mouse]",
  "16" = "Epidermal keratinocyte [Axolotl]",
  "17" = "CT progenitor [Xenopus tadpole]",
  "22" = "CT fibroblast B [Xenopus tadpole]",
  "23" = "Erythrocyte [Xenopus]",
  "24" = "CT fibroblast [Xenopus froglet]",
  "26" = "Macrophage [Mouse]",
  "27" = "Epithelial [Xenopus tadpole]",
  "28" = "Basal keratinocyte [Axolotl non-regen]"
)

# For mixed clusters: majority vote from Xenopus/axolotl cell_type column
# Collapse similar labels into shared categories
simplify_ct <- function(x) {
  x <- as.character(x)
  # Fibroblast / CT
  x <- gsub("CT_fibroblast_[A-Z]|Fibroblast$", "CT fibroblast", x)
  x <- gsub("CT_progenitor", "CT progenitor", x)
  x <- gsub("CT_cartilage_assoc|CT_perichondrial", "CT cartilage-assoc", x)
  x <- gsub("Mesenchymal Progenitor Cell", "Mesenchymal progenitor", x)
  x <- gsub("Regenerative Blastema Fibroblast", "Blastema fibroblast", x)
  # Keratinocyte / epithelial
  x <- gsub("Epidermal Keratinocyte|Epithelial Cell|Epithelial_mixed", "Keratinocyte/Epithelial", x)
  x <- gsub("Keratinocyte$", "Keratinocyte/Epithelial", x)
  x <- gsub("Basal Keratinocyte", "Basal keratinocyte", x)
  x <- gsub("Mucin-Secreting Epithelial Cell", "Mucin-secreting epithelial", x)
  # Immune
  x <- gsub("Antimicrobial Macrophage|Macrophage_2|CT_fibroblast_D", "Macrophage", x)
  x <- gsub("CD8\\+ T Cell|B Cell", "Lymphocyte", x)
  # Proliferating
  x <- gsub("Proliferating_[A-Z]", "Proliferating", x)
  # Vascular
  x <- gsub("Vascular Endothelial Cell", "Endothelial", x)
  # Neural
  x <- gsub("Schwann Cell Precursor", "Neural crest/Schwann", x)
  x <- gsub("Neural_crest", "Neural crest/Schwann", x)
  # Smooth muscle
  x <- gsub("Smooth_muscle", "Smooth muscle", x)
  x
}

message("Computing majority cell type per cluster...")
cluster_labels <- character(length(levels(atlas$seurat_clusters)))
names(cluster_labels) <- levels(atlas$seurat_clusters)

for (cl in levels(atlas$seurat_clusters)) {
  if (cl %in% names(species_specific_labels)) {
    cluster_labels[cl] <- species_specific_labels[cl]
    next
  }
  sub_ct <- atlas$cell_type[atlas$seurat_clusters == cl]
  sub_ct <- sub_ct[!is.na(sub_ct) & nzchar(as.character(sub_ct))]
  if (length(sub_ct) == 0) {
    cluster_labels[cl] <- paste0("Cluster ", cl)
    next
  }
  simple <- simplify_ct(sub_ct)
  majority <- names(which.max(table(simple)))
  cluster_labels[cl] <- majority
}

# Print the full label assignment
cat("\nCluster label assignments:\n")
for (cl in levels(atlas$seurat_clusters)) {
  n <- sum(atlas$seurat_clusters == cl)
  cat(sprintf("  Cluster %2s (n=%4d): %s\n", cl, n, cluster_labels[cl]))
}

# Add to atlas metadata — unname() prevents Seurat v5 from matching by cluster ID
atlas$cell_type_atlas <- unname(cluster_labels[as.character(atlas$seurat_clusters)])

# Save annotated atlas back to disk
message("\nSaving annotated atlas to ", ATLAS_RDS, "...")
saveRDS(atlas, ATLAS_RDS)
message("  Done.")

# =============================================================================
# 3. Annotated UMAP
# =============================================================================

pt_size <- pmin(0.5, 30000 / ncol(atlas))

p_annot <- DimPlot(atlas, reduction = "umap.harmony",
                   group.by = "cell_type_atlas",
                   label = TRUE, repel = TRUE,
                   label.size = 2.8, pt.size = pt_size * 0.7) +
  ggtitle("Cross-species atlas — cell type annotations",
          subtitle = sprintf("%d clusters | species-specific labels from marker analysis",
                             length(levels(atlas$seurat_clusters)))) +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave(file.path(FIG_DIR, "atlas_umap_annotated.png"),
       p_annot, width = 12, height = 8, dpi = 200)
message("Saved atlas_umap_annotated.png")

# =============================================================================
# 4. Load BMS bulk RNA-seq
# =============================================================================

message("\nLoading pseudobulk cache for BMS bulk data...")
cache    <- readRDS(CACHE_RDS)
bulk_mat <- cache$bulk_mat   # genes x 8 BMS/EtOH conditions
rm(cache); gc()

# Subset to BMS conditions only (exclude EtOH controls for the primary figure)
bms_cols <- grep("^BMS_|^Pretreatment$|^Reamputation$", colnames(bulk_mat), value = TRUE)
bms_mat  <- bulk_mat[, bms_cols, drop = FALSE]
cat(sprintf("BMS conditions: %s\n", paste(bms_cols, collapse = ", ")))

# =============================================================================
# 5. Per-cluster pseudobulk from atlas (mean log-normalized expression)
# =============================================================================

message("Computing per-cluster pseudobulk from atlas...")

# Genes present in both atlas and BMS bulk
shared_genes <- intersect(rownames(atlas), rownames(bulk_mat))
cat(sprintf("Shared genes (atlas ∩ BMS bulk): %d\n", length(shared_genes)))

# Get log-normalized counts (layer "data" after JoinLayers)
atlas <- JoinLayers(atlas)
expr_mat <- GetAssayData(atlas, assay = "RNA", layer = "data")
expr_mat <- expr_mat[shared_genes, ]

# Mean per cluster
clusters <- levels(atlas$seurat_clusters)
pseudo_cluster <- matrix(0, nrow = length(shared_genes), ncol = length(clusters),
                          dimnames = list(shared_genes, clusters))
for (cl in clusters) {
  cell_idx <- which(atlas$seurat_clusters == cl)
  pseudo_cluster[, cl] <- Matrix::rowMeans(expr_mat[, cell_idx, drop = FALSE])
}
cat(sprintf("Pseudobulk matrix: %d genes x %d clusters\n",
            nrow(pseudo_cluster), ncol(pseudo_cluster)))

rm(expr_mat); gc()

# =============================================================================
# 6. Spearman correlation: cluster x BMS condition
# =============================================================================

message("Computing Spearman correlations (cluster x BMS condition)...")

bms_shared <- bulk_mat[shared_genes, , drop = FALSE]

cor_mat <- matrix(NA_real_, nrow = length(clusters), ncol = ncol(bulk_mat),
                   dimnames = list(clusters, colnames(bulk_mat)))
for (cl in clusters) {
  for (cond in colnames(bulk_mat)) {
    cor_mat[cl, cond] <- cor(pseudo_cluster[, cl], bms_shared[, cond],
                              method = "spearman", use = "complete.obs")
  }
}

cat("\nSpearman rho — cluster x BMS condition (top rows):\n")
print(round(cor_mat, 3))

# Per-cluster: best BMS condition and its rho
bms_sub       <- cor_mat[, bms_cols, drop = FALSE]
which_max_col <- apply(bms_sub, 1, which.max)
best_bms_rho  <- apply(bms_sub, 1, max)
best_bms_cond <- bms_cols[which_max_col]
names(best_bms_rho)  <- clusters
names(best_bms_cond) <- clusters

cat("\nBest BMS match per cluster:\n")
for (cl in clusters) {
  cat(sprintf("  Cluster %2s  %-40s  rho=%.3f  best=%s\n",
              cl, cluster_labels[cl], best_bms_rho[cl], best_bms_cond[cl]))
}

# =============================================================================
# 7. UMAP colored by max BMS Spearman rho
# =============================================================================

message("\nGenerating BMS correlation UMAPs...")

# Add per-cell rho based on its cluster
atlas$bms_max_rho   <- unname(best_bms_rho[as.character(atlas$seurat_clusters)])
atlas$bms_best_cond <- unname(best_bms_cond[as.character(atlas$seurat_clusters)])

# Ensure bms_max_rho is numeric
atlas$bms_max_rho <- as.numeric(atlas$bms_max_rho)

# UMAP colored by max BMS rho (continuous)
# Limits set to actual data range; midpoint at mean so colors span the full gradient
rho_min <- min(atlas$bms_max_rho)
rho_max <- max(atlas$bms_max_rho)
rho_mid <- mean(c(rho_min, rho_max))

p_rho <- FeaturePlot(atlas, reduction = "umap.harmony",
                      features = "bms_max_rho", pt.size = pt_size) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#B2182B",
    midpoint = rho_mid,
    limits   = c(rho_min, rho_max),
    name     = "Spearman rho\n(best BMS match)"
  ) +
  ggtitle("BMS bulk similarity — max Spearman rho per cluster",
          subtitle = "Higher rho = more similar to BMS bulk RNA-seq") +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave(file.path(FIG_DIR, "atlas_umap_bms_rho.png"),
       p_rho, width = 9, height = 6, dpi = 200)
message("Saved atlas_umap_bms_rho.png")

# UMAP colored by best BMS condition (categorical)
bms_cond_pal <- c(
  BMS_acute      = "#B71C1C",
  BMS_24hrs      = "#C62828",
  Pretreatment   = "#E53935",
  Reamputation   = "#EF5350",
  EtOH_acute     = "#1565C0",
  EtOH_24hrs     = "#1976D2",
  EtOH_pretreat  = "#1E88E5",
  EtOH_reamp     = "#42A5F5"
)

p_cond <- DimPlot(atlas, reduction = "umap.harmony",
                   group.by = "bms_best_cond",
                   cols = bms_cond_pal,
                   pt.size = pt_size, shuffle = FALSE) +
  ggtitle("Best-matching BMS condition per cluster",
          subtitle = "Each cluster colored by the BMS/EtOH condition with highest Spearman rho") +
  theme(plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave(file.path(FIG_DIR, "atlas_umap_bms_condition.png"),
       p_cond, width = 10, height = 6, dpi = 200)
message("Saved atlas_umap_bms_condition.png")

# =============================================================================
# 8. Heatmap: cluster x all BMS+EtOH conditions
# =============================================================================

message("Generating cluster x BMS correlation heatmap...")

# Row labels: cluster number + short cell type
row_labels <- paste0(clusters, ": ", cluster_labels[clusters])

# Order clusters by hierarchical clustering of BMS conditions
hc_order <- hclust(dist(cor_mat))$order

cor_long <- data.frame(
  cluster   = factor(rep(row_labels, ncol(cor_mat)),
                     levels = row_labels[hc_order]),
  condition = rep(colnames(cor_mat), each = nrow(cor_mat)),
  rho       = as.vector(cor_mat)
)

# Condition order: BMS first, then EtOH
cond_order <- c(bms_cols,
                grep("^EtOH_", colnames(bulk_mat), value = TRUE))
cor_long$condition <- factor(cor_long$condition, levels = cond_order)

p_heat <- ggplot(cor_long, aes(x = condition, y = cluster, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.2, color = "black") +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#B2182B",
    midpoint = mean(range(cor_mat)),
    limits   = range(cor_mat),
    name     = "Spearman rho"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title    = "Atlas cluster x BMS bulk Spearman correlation",
       subtitle  = sprintf("%d shared genes | per-cluster pseudobulk vs bulk RNA-seq",
                           length(shared_genes)),
       x = "BMS / EtOH condition",
       y = "Atlas cluster") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x  = element_text(angle = 40, hjust = 1, size = 9),
        axis.text.y  = element_text(size  = 7),
        plot.title   = element_text(size  = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        panel.grid   = element_blank())

ggsave(file.path(FIG_DIR, "atlas_bms_correlation_heatmap.png"),
       p_heat, width = 10, height = 11, dpi = 200)
message("Saved atlas_bms_correlation_heatmap.png")

# =============================================================================
# 9. Combined panel: annotated UMAP + BMS rho UMAP side-by-side
# =============================================================================

p_combined <- p_annot + p_rho +
  plot_annotation(title    = "Cross-species atlas — cell types and BMS similarity",
                  subtitle = "Left: annotated clusters | Right: Spearman rho to best BMS bulk condition")

ggsave(file.path(FIG_DIR, "atlas_annotated_bms_combined.png"),
       p_combined, width = 22, height = 8, dpi = 200)
message("Saved atlas_annotated_bms_combined.png")

message("\nDone. Figures in ", FIG_DIR)
