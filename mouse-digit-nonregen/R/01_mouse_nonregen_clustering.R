# =============================================================================
# 01 — Mouse digit (non-regenerative): QC, merge, normalise, cluster, UMAP
#
# Input:  data/raw/  (11 samples — run 00_download_data.R first)
# Output: results/mouse_nonregen_merged.rds
#
# Storer et al. 2020 (Dev Cell) — GSE135985  PMID 31902657
# Same Lehoczky lab as Johnson, Masias & Lehoczky 2020/GSE143888 (processed in mouse-digit/).
#
# Two conditions:
#   regenerative    : distal digit tip amputations (0, 7, 10, 14, 28, 56 dpa)
#   non-regenerative: proximal amputations — fail to regrow (10, 14 dpa)
#
# Metadata columns:
#   condition : "regenerative" | "non-regenerative"
#   timepoint : "uninjured","7dpa","10dpa","14dpa","28dpa","56dpa" (regen)
#               "10dpa","14dpa" (non-regen — stored with condition prefix in sample)
#   sample    : individual library name (used by Harmony)
#   species   : "Mouse"
# =============================================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# Sample name -> (condition, timepoint) mapping
SAMPLE_META <- list(
  Regen_7dpa         = list(cond = "regenerative",     tp = "7dpa"),
  Regen_10dpa        = list(cond = "regenerative",     tp = "10dpa"),
  Regen_14dpa_A      = list(cond = "regenerative",     tp = "14dpa"),
  Regen_14dpa_B      = list(cond = "regenerative",     tp = "14dpa"),
  Regen_28dpa        = list(cond = "regenerative",     tp = "28dpa"),
  Regen_56dpa        = list(cond = "regenerative",     tp = "56dpa"),
  Regen_uninjured_A  = list(cond = "regenerative",     tp = "uninjured"),
  Regen_uninjured_B  = list(cond = "regenerative",     tp = "uninjured"),
  NonRegen_10dpa_A   = list(cond = "non-regenerative", tp = "10dpa"),
  NonRegen_10dpa_B   = list(cond = "non-regenerative", tp = "10dpa"),
  NonRegen_14dpa     = list(cond = "non-regenerative", tp = "14dpa")
)

COND_COLORS <- c("regenerative" = "#1565C0", "non-regenerative" = "#B71C1C")
TP_ORDER    <- c("uninjured","7dpa","10dpa","14dpa","28dpa","56dpa")

# =============================================================================
# 1. LOAD RAW 10x DATA
# =============================================================================

load_10x <- function(path, project_name) {
  counts <- Read10X(data.dir = path)
  if (is.list(counts)) counts <- counts[[1]]
  CreateSeuratObject(counts = counts, project = project_name, min.cells = 3)
}

raw <- file.path("data", "raw")
message("Loading GSE135985 samples...")

objs <- lapply(names(SAMPLE_META), function(nm) {
  path <- file.path(raw, nm)
  if (!dir.exists(path)) {
    message("  MISSING: ", nm, " — run 00_download_data.R first")
    return(NULL)
  }
  message("  Loading ", nm)
  load_10x(path, nm)
})
names(objs) <- names(SAMPLE_META)
objs <- Filter(Negate(is.null), objs)

message("Cells per sample: ")
for (nm in names(objs)) message("  ", nm, ": ", ncol(objs[[nm]]))

# =============================================================================
# 2. PER-SAMPLE QC
# =============================================================================

add_mito <- function(seu) {
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

objs <- lapply(objs, add_mito)

p_qc <- VlnPlot(
  merge(objs[[1]], y = objs[-1]),
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  ncol = 3, pt.size = 0
)
pdf(file.path(results_dir, "mouse_nonregen_qc_violin_prefilter.pdf"), width = 22, height = 4)
print(p_qc)
dev.off()
message("QC violin saved — inspect before adjusting thresholds if needed")

objs <- lapply(objs, function(seu) {
  subset(seu, nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 20)
})
message("Post-QC cells: ", paste(sapply(objs, ncol), collapse = " | "))

# =============================================================================
# 3. ADD METADATA
# =============================================================================

objs <- lapply(names(objs), function(nm) {
  seu <- objs[[nm]]
  m   <- SAMPLE_META[[nm]]
  seu$condition <- m$cond
  seu$timepoint <- m$tp
  seu$sample    <- nm
  seu$species   <- "Mouse"
  seu
})
names(objs) <- names(SAMPLE_META)[seq_along(objs)]

# =============================================================================
# 4. MERGE + HARMONY INTEGRATION
# =============================================================================

message("Merging all samples...")
merged <- merge(
  objs[[1]],
  y           = objs[-1],
  add.cell.ids = names(objs),
  project     = "MouseDigitRegenVsNonRegen"
)

rm(objs); gc()

merged$timepoint <- factor(merged$timepoint, levels = TP_ORDER)

message("Total merged cells: ", ncol(merged))
print(table(merged$condition, merged$timepoint))

# =============================================================================
# 5. NORMALISE, HVG, SCALE, PCA
# =============================================================================

message("Normalising and running PCA...")
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 2000)
merged <- ScaleData(merged, features = VariableFeatures(merged))
merged <- RunPCA(merged, npcs = 30)
gc()

pdf(file.path(results_dir, "mouse_nonregen_elbow.pdf"), width = 6, height = 4)
print(ElbowPlot(merged, ndims = 30))
dev.off()

message("Running Harmony (batch correction by sample)...")
merged <- RunHarmony(
  merged,
  group.by.vars  = "sample",
  reduction      = "pca",
  reduction.save = "harmony",
  verbose        = FALSE
)

# =============================================================================
# 6. CLUSTERING + UMAP
# =============================================================================

merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:20)
merged <- FindClusters(merged, resolution = 0.3)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:20,
                   reduction.name = "umap_harmony",
                   min.dist = 0.3, spread = 1,
                   n.neighbors = 30, n.epochs = 200)

merged@misc$default_reduction <- "umap_harmony"

message("Clusters: ", length(unique(Idents(merged))))
message("Total cells: ", ncol(merged))

pdf(file.path(results_dir, "mouse_nonregen_umap_clusters.pdf"), width = 8, height = 6)
print(DimPlot(merged, reduction = "umap_harmony", label = TRUE, repel = TRUE) +
        ggtitle("Mouse digit (regen + non-regen) — Seurat clusters (res=0.3)"))
dev.off()

pdf(file.path(results_dir, "mouse_nonregen_umap_condition.pdf"), width = 10, height = 6)
print(DimPlot(merged, reduction = "umap_harmony", group.by = "condition",
              cols = COND_COLORS, pt.size = 0.4) +
        ggtitle("Mouse digit — regenerative vs non-regenerative"))
dev.off()

pdf(file.path(results_dir, "mouse_nonregen_umap_timepoint.pdf"), width = 12, height = 6)
print(DimPlot(merged, reduction = "umap_harmony", group.by = "timepoint",
              pt.size = 0.4) +
        ggtitle("Mouse digit — by timepoint"))
dev.off()

# =============================================================================
# 7. SAVE (strip heavy layers)
# =============================================================================

merged[["RNA"]]$scale.data <- NULL
for (l in Layers(merged)[startsWith(Layers(merged), "counts.")])
  LayerData(merged, layer = l) <- NULL
gc()

saveRDS(merged, file.path(results_dir, "mouse_nonregen_merged.rds"))
message("Saved results/mouse_nonregen_merged.rds")
message("Run 02_mouse_nonregen_annotation.R next — inspect UMAPs and QC first.")
