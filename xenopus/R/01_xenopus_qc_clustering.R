# =============================================================================
# 01 — Xenopus laevis: QC, merge, normalise, cluster, UMAP
# Input:  data/raw/Xen_BL0dpa/, Xen_BL3dpa/, Xen_Pool_BL7_10_14dpa/,
#                  Xen_BL14dpa/, Xen_Pool_BL14_20_52dpa/
# Output: results/xen_merged.rds
#
# Five timepoints spanning 0-52 dpa:
#   0dpa     — intact/pre-regen baseline
#   3dpa     — early wound response
#   7-14dpa  — early blastema (pooled)
#   14dpa    — mid-blastema (standalone replicate)
#   14-52dpa — late blastema / maturation (pooled)
#
# The 14 dpa standalone and 14-20-52 dpa pool extend the original three
# samples to capture late-stage GOI (Hedgehog pathway, patterning genes)
# that are expected to activate at blastema specification stages.
# =============================================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

TIMEPOINT_ORDER <- c("0dpa","3dpa","7-14dpa","14dpa","14-52dpa")

# =============================================================================
# 1. LOAD RAW 10x DATA
# =============================================================================

load_10x <- function(path, project_name) {
  counts <- Read10X(data.dir = path)
  CreateSeuratObject(counts = counts, project = project_name, min.cells = 3)
}

message("Loading Xenopus 10x data...")
xen_0dpa    <- load_10x(file.path("data","raw","Xen_BL0dpa"),              "Xen_0dpa")
xen_3dpa    <- load_10x(file.path("data","raw","Xen_BL3dpa"),              "Xen_3dpa")
xen_pool    <- load_10x(file.path("data","raw","Xen_Pool_BL7_10_14dpa"),   "Xen_pool")
xen_14dpa   <- load_10x(file.path("data","raw","Xen_BL14dpa"),             "Xen_14dpa")
xen_late    <- load_10x(file.path("data","raw","Xen_Pool_BL14_20_52dpa"),  "Xen_late")

message(sprintf("Cells loaded — 0dpa: %d | 3dpa: %d | pool: %d | 14dpa: %d | late: %d",
                ncol(xen_0dpa), ncol(xen_3dpa), ncol(xen_pool),
                ncol(xen_14dpa), ncol(xen_late)))

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

xen_0dpa  <- add_mito(xen_0dpa)
xen_3dpa  <- add_mito(xen_3dpa)
xen_pool  <- add_mito(xen_pool)
xen_14dpa <- add_mito(xen_14dpa)
xen_late  <- add_mito(xen_late)

# QC plots before filtering — inspect before setting thresholds
p_qc <- VlnPlot(
  merge(xen_0dpa, list(xen_3dpa, xen_pool, xen_14dpa, xen_late)),
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  ncol = 3, pt.size = 0
)
pdf(file.path(results_dir, "xenopus_qc_violin_prefilter.pdf"), width = 18, height = 4)
print(p_qc)
dev.off()
message("QC violin saved — inspect before adjusting thresholds below")

# Thresholds for original three samples established from violin plots.
# Thresholds for xen_14dpa and xen_late are conservative starting points —
# inspect xenopus_qc_violin_prefilter.pdf and adjust if needed.
xen_0dpa  <- subset(xen_0dpa,  nFeature_RNA > 500 & nFeature_RNA < 8000  & percent.mt < 20)
xen_3dpa  <- subset(xen_3dpa,  nFeature_RNA > 500 & nFeature_RNA < 12000 & percent.mt < 20)
xen_pool  <- subset(xen_pool,  nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)
xen_14dpa <- subset(xen_14dpa, nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)
xen_late  <- subset(xen_late,  nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)

message(sprintf(
  "Post-QC cells — 0dpa: %d | 3dpa: %d | 7-14dpa: %d | 14dpa: %d | 14-52dpa: %d",
  ncol(xen_0dpa), ncol(xen_3dpa), ncol(xen_pool),
  ncol(xen_14dpa), ncol(xen_late)
))
# Original expected: 3,869 | 5,251 | 10,828

# =============================================================================
# 3. ADD METADATA
# =============================================================================

xen_0dpa$timepoint  <- "0dpa";     xen_0dpa$sample  <- "Xen_BL0dpa"
xen_3dpa$timepoint  <- "3dpa";     xen_3dpa$sample  <- "Xen_BL3dpa"
xen_pool$timepoint  <- "7-14dpa";  xen_pool$sample  <- "Xen_Pool_7-14dpa"
xen_14dpa$timepoint <- "14dpa";    xen_14dpa$sample <- "Xen_BL14dpa"
xen_late$timepoint  <- "14-52dpa"; xen_late$sample  <- "Xen_Pool_14-52dpa"

# Direct assignment required — R for-loop copies don't write back to originals
xen_0dpa$species  <- "Xenopus"
xen_3dpa$species  <- "Xenopus"
xen_pool$species  <- "Xenopus"
xen_14dpa$species <- "Xenopus"
xen_late$species  <- "Xenopus"

# =============================================================================
# 4. MERGE + HARMONY INTEGRATION
# The original 3 samples (GSM5045xxx) and the 2 new samples (GSM5057xxx) are
# from different GEO submissions and require batch correction. Harmony is run
# on sample identity (group.by.vars = "sample"), matching the approach used
# for the Leigh 2018 axolotl multi-study integration.
# =============================================================================

message("Merging Xenopus objects...")
xen_merged <- merge(
  xen_0dpa,
  y            = list(xen_3dpa, xen_pool, xen_14dpa, xen_late),
  add.cell.ids = c("0dpa","3dpa","pool","14dpa","late"),
  project      = "XenopusRegen"
)

xen_merged <- JoinLayers(xen_merged)
xen_merged$timepoint <- factor(xen_merged$timepoint, levels = TIMEPOINT_ORDER)

message("Total merged cells: ", ncol(xen_merged))

# =============================================================================
# 5. NORMALISE, HVG, SCALE, PCA
# =============================================================================

message("Normalising, integrating, and clustering...")
xen_merged <- NormalizeData(xen_merged)
xen_merged <- FindVariableFeatures(xen_merged, nfeatures = 2000)
xen_merged <- ScaleData(xen_merged)
xen_merged <- RunPCA(xen_merged, npcs = 30)

pdf(file.path(results_dir, "xenopus_elbow_plot.pdf"), width = 6, height = 4)
print(ElbowPlot(xen_merged, ndims = 30))
dev.off()

message("Running Harmony (batch correction by sample)...")
xen_merged <- RunHarmony(
  xen_merged,
  group.by.vars  = "sample",
  reduction      = "pca",
  reduction.save = "harmony",
  verbose        = FALSE
)

# =============================================================================
# 6. CLUSTERING + UMAP on Harmony embedding
# =============================================================================

xen_merged <- FindNeighbors(xen_merged, reduction = "harmony", dims = 1:20)
xen_merged <- FindClusters(xen_merged, resolution = 0.3)
xen_merged <- RunUMAP(xen_merged, reduction = "harmony", dims = 1:20,
                       reduction.name = "umap_harmony",
                       min.dist = 0.3, spread = 1,
                       n.neighbors = 30, n.epochs = 200)

xen_merged@misc$default_reduction <- "umap_harmony"

message("Clusters: ", length(unique(Idents(xen_merged))))
message("Total cells: ", ncol(xen_merged))
print(table(xen_merged$timepoint))

pdf(file.path(results_dir, "xenopus_umap_clusters.pdf"), width = 8, height = 6)
print(DimPlot(xen_merged, reduction = "umap_harmony", label = TRUE, repel = TRUE) +
        ggtitle("Xenopus - Seurat clusters (res=0.3)"))
dev.off()

pdf(file.path(results_dir, "xenopus_umap_timepoint.pdf"), width = 10, height = 6)
print(DimPlot(xen_merged, reduction = "umap_harmony", group.by = "timepoint", pt.size = 0.4,
              cols = c("0dpa"     = "#2E7D32",
                       "3dpa"     = "#00838F",
                       "7-14dpa"  = "#1565C0",
                       "14dpa"    = "#6A1B9A",
                       "14-52dpa" = "#B71C1C")) +
        ggtitle("Xenopus - by timepoint"))
dev.off()

saveRDS(xen_merged, file.path(results_dir, "xen_merged.rds"))
message("Saved results/xen_merged.rds")
