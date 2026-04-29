# =============================================================================
# 01 — Mouse digit tip: QC, merge, normalise, cluster, UMAP
#
# Input:  data/raw/  (5 samples — run 00_download_data.R first)
# Output: results/mouse_merged.rds
#
# Storer et al. 2020 (Dev Cell) — GSE143888
# Adult mouse digit tip regeneration (distal phalanx amputation)
#
# Samples:
#   Mouse_0dpa   — unamputated control      (GSM4276219)
#   Mouse_11dpa  — 11 days post-amputation  (GSM4276220)
#   Mouse_12dpa  — 12 days post-amputation  (GSM4276221)
#   Mouse_14dpa  — 14 days post-amputation  (GSM4276222)
#   Mouse_17dpa  — 17 days post-amputation  (GSM4276223)
#
# Metadata columns:
#   timepoint : "0dpa","11dpa","12dpa","14dpa","17dpa"
#   sample    : individual library name (used by Harmony)
#   species   : "Mouse"
#   condition : "regenerative" (all samples — mouse digit tip is regeneration-competent)
# =============================================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

TIMEPOINT_ORDER <- c("0dpa","11dpa","12dpa","14dpa","17dpa")

TIMEPOINT_COLORS <- c(
  "0dpa"   = "#1B5E20",
  "11dpa"  = "#2E7D32",
  "12dpa"  = "#388E3C",
  "14dpa"  = "#66BB6A",
  "17dpa"  = "#A5D6A7"
)

# =============================================================================
# 1. LOAD RAW 10x DATA
# =============================================================================

load_10x <- function(path, project_name) {
  counts <- Read10X(data.dir = path)
  # Read10X may return a list for multi-modal data — take the first element
  if (is.list(counts)) counts <- counts[[1]]
  CreateSeuratObject(counts = counts, project = project_name, min.cells = 3)
}

raw <- file.path("data", "raw")
message("Loading Mouse digit tip 10x data...")

mouse_0dpa  <- load_10x(file.path(raw, "Mouse_0dpa"),  "Mouse_0dpa")
mouse_11dpa <- load_10x(file.path(raw, "Mouse_11dpa"), "Mouse_11dpa")
mouse_12dpa <- load_10x(file.path(raw, "Mouse_12dpa"), "Mouse_12dpa")
mouse_14dpa <- load_10x(file.path(raw, "Mouse_14dpa"), "Mouse_14dpa")
mouse_17dpa <- load_10x(file.path(raw, "Mouse_17dpa"), "Mouse_17dpa")

all_objs <- list(mouse_0dpa, mouse_11dpa, mouse_12dpa, mouse_14dpa, mouse_17dpa)
message(sprintf("Raw cells per sample: %s",
        paste(sapply(all_objs, ncol), collapse = " | ")))

# Check gene ID format — mouse 10x is usually gene symbols (Col1a1) or Ensembl
sample_genes <- head(rownames(mouse_0dpa), 10)
message("Sample gene IDs: ", paste(sample_genes, collapse = ", "))

# =============================================================================
# 2. PER-SAMPLE QC
# =============================================================================

add_mito <- function(seu) {
  # Mouse mito genes start with "mt-" (lowercase)
  mt_genes <- grep("^mt-", rownames(seu), value = TRUE, ignore.case = FALSE)
  if (length(mt_genes) == 0)
    mt_genes <- grep("^MT-", rownames(seu), value = TRUE, ignore.case = FALSE)
  if (length(mt_genes) > 0) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mt_genes)
  } else {
    seu[["percent.mt"]] <- 0
    message("  No mito genes found for ", seu@project.name, " — percent.mt set to 0")
  }
  seu
}

all_objs <- lapply(all_objs, add_mito)
list2env(setNames(all_objs,
  c("mouse_0dpa","mouse_11dpa","mouse_12dpa","mouse_14dpa","mouse_17dpa")),
  envir = .GlobalEnv)

p_qc <- VlnPlot(
  merge(mouse_0dpa, list(mouse_11dpa, mouse_12dpa, mouse_14dpa, mouse_17dpa)),
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  ncol = 3, pt.size = 0
)
pdf(file.path(results_dir, "mouse_qc_violin_prefilter.pdf"), width = 18, height = 4)
print(p_qc)
dev.off()
message("QC violin saved — inspect mouse_qc_violin_prefilter.pdf before adjusting thresholds")

# Conservative starting thresholds — inspect violin and adjust
mouse_0dpa  <- subset(mouse_0dpa,  nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
mouse_11dpa <- subset(mouse_11dpa, nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
mouse_12dpa <- subset(mouse_12dpa, nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
mouse_14dpa <- subset(mouse_14dpa, nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
mouse_17dpa <- subset(mouse_17dpa, nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)

message("Post-QC cells:")
message(sprintf("  0dpa=%d | 11dpa=%d | 12dpa=%d | 14dpa=%d | 17dpa=%d",
  ncol(mouse_0dpa), ncol(mouse_11dpa), ncol(mouse_12dpa),
  ncol(mouse_14dpa), ncol(mouse_17dpa)))

# =============================================================================
# 3. ADD METADATA
# =============================================================================

set_meta <- function(seu, timepoint, sample) {
  seu$timepoint <- timepoint
  seu$sample    <- sample
  seu$species   <- "Mouse"
  seu$condition <- "regenerative"
  seu
}

mouse_0dpa  <- set_meta(mouse_0dpa,  "0dpa",  "Mouse_0dpa")
mouse_11dpa <- set_meta(mouse_11dpa, "11dpa", "Mouse_11dpa")
mouse_12dpa <- set_meta(mouse_12dpa, "12dpa", "Mouse_12dpa")
mouse_14dpa <- set_meta(mouse_14dpa, "14dpa", "Mouse_14dpa")
mouse_17dpa <- set_meta(mouse_17dpa, "17dpa", "Mouse_17dpa")

# =============================================================================
# 4. MERGE + HARMONY INTEGRATION
# =============================================================================

message("Merging all Mouse objects...")
mouse_merged <- merge(
  mouse_0dpa,
  y = list(mouse_11dpa, mouse_12dpa, mouse_14dpa, mouse_17dpa),
  add.cell.ids = c("0dpa","11dpa","12dpa","14dpa","17dpa"),
  project = "MouseDigitRegen"
)

rm(mouse_0dpa, mouse_11dpa, mouse_12dpa, mouse_14dpa, mouse_17dpa, all_objs)
gc()

mouse_merged$timepoint <- factor(mouse_merged$timepoint, levels = TIMEPOINT_ORDER)

message("Total merged cells: ", ncol(mouse_merged))
print(table(mouse_merged$timepoint))

# =============================================================================
# 5. NORMALISE, HVG, SCALE, PCA
# =============================================================================

message("Normalising and running PCA...")
mouse_merged <- NormalizeData(mouse_merged)
mouse_merged <- FindVariableFeatures(mouse_merged, nfeatures = 2000)
mouse_merged <- ScaleData(mouse_merged, features = VariableFeatures(mouse_merged))
mouse_merged <- RunPCA(mouse_merged, npcs = 30)
gc()

pdf(file.path(results_dir, "mouse_elbow_plot.pdf"), width = 6, height = 4)
print(ElbowPlot(mouse_merged, ndims = 30))
dev.off()

message("Running Harmony (batch correction by sample)...")
mouse_merged <- RunHarmony(
  mouse_merged,
  group.by.vars  = "sample",
  reduction      = "pca",
  reduction.save = "harmony",
  verbose        = FALSE
)

# =============================================================================
# 6. CLUSTERING + UMAP
# =============================================================================

mouse_merged <- FindNeighbors(mouse_merged, reduction = "harmony", dims = 1:20)
mouse_merged <- FindClusters(mouse_merged, resolution = 0.3)
mouse_merged <- RunUMAP(mouse_merged, reduction = "harmony", dims = 1:20,
                         reduction.name = "umap_harmony",
                         min.dist = 0.3, spread = 1,
                         n.neighbors = 30, n.epochs = 200)

mouse_merged@misc$default_reduction <- "umap_harmony"

message("Clusters: ", length(unique(Idents(mouse_merged))))
message("Total cells: ", ncol(mouse_merged))

pdf(file.path(results_dir, "mouse_umap_clusters.pdf"), width = 8, height = 6)
print(DimPlot(mouse_merged, reduction = "umap_harmony", label = TRUE, repel = TRUE) +
        ggtitle("Mouse digit tip - Seurat clusters (res=0.3)"))
dev.off()

pdf(file.path(results_dir, "mouse_umap_timepoint.pdf"), width = 10, height = 6)
print(DimPlot(mouse_merged, reduction = "umap_harmony", group.by = "timepoint",
              cols = TIMEPOINT_COLORS, pt.size = 0.4) +
        ggtitle("Mouse digit tip - by timepoint"))
dev.off()

# Drop heavy layers before saving
mouse_merged[["RNA"]]$scale.data <- NULL
count_layers <- Layers(mouse_merged)[startsWith(Layers(mouse_merged), "counts.")]
for (l in count_layers) LayerData(mouse_merged, layer = l) <- NULL
gc()

saveRDS(mouse_merged, file.path(results_dir, "mouse_merged.rds"))
message("Saved results/mouse_merged.rds")
message("Run 02_mouse_annotation.R next — inspect mouse_qc_violin_prefilter.pdf first.")
