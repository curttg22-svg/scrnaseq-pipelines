# =============================================================================
# 01 — Xenopus laevis: QC, merge, normalise, cluster, UMAP
# Input:  data/raw/Xen_BL0dpa/, Xen_BL3dpa/, Xen_Pool_BL7_10_14dpa/
# Output: results/xen_merged.rds
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# =============================================================================
# 1. LOAD RAW 10x DATA
# =============================================================================

load_10x <- function(path, project_name) {
  counts <- Read10X(data.dir = path)
  CreateSeuratObject(counts = counts, project = project_name, min.cells = 3)
}

message("Loading Xenopus 10x data...")
xen_0dpa  <- load_10x(file.path("data","raw","Xen_BL0dpa"),          "Xen_0dpa")
xen_3dpa  <- load_10x(file.path("data","raw","Xen_BL3dpa"),          "Xen_3dpa")
xen_pool  <- load_10x(file.path("data","raw","Xen_Pool_BL7_10_14dpa"),"Xen_pool")

# =============================================================================
# 2. PER-SAMPLE QC
# =============================================================================

add_mito <- function(seu) {
  # Xenopus mito genes are ALL-CAPS (ND1-6, CYTB, COX1-3, ATP6/8)
  # Try uppercase first, fall back to lowercase, set 0 if none found
  mt_genes <- grep("^MT-|^ND[0-9]|^CYTB|^COX[0-9]|^ATP[68]",
                   rownames(seu), value = TRUE, ignore.case = FALSE)
  if (length(mt_genes) == 0)
    mt_genes <- grep("^mt-", rownames(seu), value = TRUE, ignore.case = FALSE)
  if (length(mt_genes) > 0) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mt_genes)
  } else {
    seu[["percent.mt"]] <- 0
    message("  No mito genes found for ", seu@project.name, " — percent.mt set to 0")
  }
  seu
}

xen_0dpa <- add_mito(xen_0dpa)
xen_3dpa <- add_mito(xen_3dpa)
xen_pool <- add_mito(xen_pool)

# QC plots before filtering
p_qc <- VlnPlot(merge(xen_0dpa, list(xen_3dpa, xen_pool)),
                 features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                 ncol = 3, pt.size = 0)
pdf(file.path(results_dir, "xenopus_qc_violin_prefilter.pdf"), width = 12, height = 4)
print(p_qc)
dev.off()

# Thresholds established by inspecting violin plots
xen_0dpa <- subset(xen_0dpa, nFeature_RNA > 500 & nFeature_RNA < 8000  & percent.mt < 20)
xen_3dpa <- subset(xen_3dpa, nFeature_RNA > 500 & nFeature_RNA < 12000 & percent.mt < 20)
xen_pool <- subset(xen_pool, nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)

message(sprintf("Post-QC cells — 0dpa: %d | 3dpa: %d | pool: %d",
                ncol(xen_0dpa), ncol(xen_3dpa), ncol(xen_pool)))
# Expected: 3,869 | 5,251 | 10,828

# =============================================================================
# 3. ADD METADATA
# =============================================================================

xen_0dpa$timepoint <- "0dpa";     xen_0dpa$sample <- "Xen_BL0dpa"
xen_3dpa$timepoint <- "3dpa";     xen_3dpa$sample <- "Xen_BL3dpa"
xen_pool$timepoint <- "7-14dpa";  xen_pool$sample <- "Xen_Pool_7-14dpa"

# Direct assignment required — R for-loop copies don't write back to originals
xen_0dpa$species <- "Xenopus"
xen_3dpa$species <- "Xenopus"
xen_pool$species <- "Xenopus"

# =============================================================================
# 4. MERGE — no Harmony needed; timepoint separation is biological
# =============================================================================

message("Merging Xenopus objects...")
xen_merged <- merge(xen_0dpa,
                     y         = list(xen_3dpa, xen_pool),
                     add.cell.ids = c("0dpa","3dpa","pool"),
                     project   = "XenopusRegen")

xen_merged <- JoinLayers(xen_merged)   # required in Seurat v5 after merge

# =============================================================================
# 5. NORMALISE, HVG, SCALE, PCA
# =============================================================================

message("Normalising and clustering...")
xen_merged <- NormalizeData(xen_merged)
xen_merged <- FindVariableFeatures(xen_merged, nfeatures = 2000)
xen_merged <- ScaleData(xen_merged)
xen_merged <- RunPCA(xen_merged, npcs = 30)

# ElbowPlot -> PC 20 captures most variance
pdf(file.path(results_dir, "xenopus_elbow_plot.pdf"), width = 6, height = 4)
print(ElbowPlot(xen_merged, ndims = 30))
dev.off()

# =============================================================================
# 6. CLUSTERING + UMAP
# =============================================================================

xen_merged <- FindNeighbors(xen_merged, dims = 1:20)
xen_merged <- FindClusters(xen_merged,  resolution = 0.3)  # -> 20 clusters
xen_merged <- RunUMAP(xen_merged, dims = 1:20,
                       min.dist = 0.3, spread = 1,
                       n.neighbors = 30, n.epochs = 200)

message("Clusters: ", length(unique(Idents(xen_merged))))
message("Total cells: ", ncol(xen_merged))

pdf(file.path(results_dir, "xenopus_umap_clusters.pdf"), width = 8, height = 6)
print(DimPlot(xen_merged, label = TRUE, repel = TRUE) +
        ggtitle("Xenopus — Seurat clusters (res=0.3)"))
dev.off()

pdf(file.path(results_dir, "xenopus_umap_timepoint.pdf"), width = 8, height = 6)
print(DimPlot(xen_merged, group.by = "timepoint", pt.size = 0.4) +
        ggtitle("Xenopus — by timepoint"))
dev.off()

saveRDS(xen_merged, file.path(results_dir, "xen_merged.rds"))
message("Saved results/xen_merged.rds")
