# =============================================================================
# 03 — Axolotl (A. mexicanum): QC, normalise, cluster, UMAP
# Input:  data/raw/Axo_BL11dpa_rep1/
# Output: results/axo_final.rds
#
# Note: dataset is CT fate-map enriched — ~97% of cells are mesenchymal/CT
#       lineage by experimental design. This is expected, not a limitation.
# Note: axolotl rep2 (GSM5045049) excluded — only 1,765 cells vs 5,391 in
#       rep1, with a batch effect not correctable by Harmony.
# =============================================================================

library(Seurat)
library(ggplot2)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# =============================================================================
# 1. LOAD + QC
# =============================================================================

message("Loading axolotl 10x data...")
axo_counts <- Read10X(data.dir = file.path("data","raw","Axo_BL11dpa_rep1"))
axo_final  <- CreateSeuratObject(axo_counts, project = "AxolotlRegen", min.cells = 3)

# Mitochondrial genes — axolotl mt genes don't use "mt-" prefix reliably;
# use the AMEX IDs for mt-encoded genes if available, else skip
# (percent.mt used only for outlier cell detection)
mt_genes <- grep("^mt-|^MT-", rownames(axo_final), value = TRUE, ignore.case = TRUE)
if (length(mt_genes) > 0) {
  axo_final[["percent.mt"]] <- PercentageFeatureSet(axo_final, features = mt_genes)
} else {
  axo_final[["percent.mt"]] <- 0
  message("No mt- genes found in AMEX reference — percent.mt set to 0")
}

# QC violin before filter
pdf(file.path(results_dir, "axolotl_qc_violin_prefilter.pdf"), width = 9, height = 4)
print(VlnPlot(axo_final,
              features = c("nFeature_RNA","nCount_RNA","percent.mt"),
              ncol = 3, pt.size = 0))
dev.off()

# 93–97% of barcodes are empty droplets — aggressive lower threshold required
axo_final <- subset(axo_final,
                     nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 20)

message("Post-QC cells: ", ncol(axo_final))  # expected ~5,391

axo_final$species   <- "Axolotl"
axo_final$timepoint <- "11dpa"
axo_final$sample    <- "Axo_BL11dpa_rep1"

# =============================================================================
# 2. NORMALISE, HVG, SCALE, PCA
# =============================================================================

message("Normalising and clustering axolotl data...")
axo_final <- JoinLayers(axo_final)
axo_final <- NormalizeData(axo_final)
axo_final <- FindVariableFeatures(axo_final, nfeatures = 3000)
axo_final <- ScaleData(axo_final)
axo_final <- RunPCA(axo_final, npcs = 30)

pdf(file.path(results_dir, "axolotl_elbow_plot.pdf"), width = 6, height = 4)
print(ElbowPlot(axo_final, ndims = 30))
dev.off()

# =============================================================================
# 3. CLUSTERING + UMAP
# ElbowPlot -> PC15
# =============================================================================

axo_final <- FindNeighbors(axo_final, dims = 1:15)
axo_final <- FindClusters(axo_final,  resolution = 0.5)  # -> 7 clusters
axo_final <- RunUMAP(axo_final, dims = 1:15, min.dist = 0.3, spread = 1)

message("Clusters: ", length(unique(Idents(axo_final))))

pdf(file.path(results_dir, "axolotl_umap_clusters.pdf"), width = 7, height = 6)
print(DimPlot(axo_final, label = TRUE, repel = TRUE) +
        ggtitle("Axolotl — Seurat clusters (res=0.5)"))
dev.off()

saveRDS(axo_final, file.path(results_dir, "axo_final.rds"))
message("Saved results/axo_final.rds")
