# =============================================================================
# 01 — Xenopus laevis: QC, merge, normalise, cluster, UMAP
#
# Input:  data/raw/  (10 samples — run 00_download_data.R first)
# Output: results/xen_merged.rds
#
# Two experimental groups from Lin et al. 2021 (GSE165901):
#
#   BL  (blastema, NON-regenerative froglet post-amputation):
#     Xen_BL0dpa              — 0 dpa     (GSM5045045)
#     Xen_FACS_BL0dpa         — 0 dpa     (GSM5057655, FACS-sorted, same-animal supplement)
#     Xen_BL3dpa              — 3 dpa     (GSM5045046)
#     Xen_FACS_BL3dpa         — 3 dpa     (GSM5057656, FACS-sorted, same-animal supplement)
#     Xen_Pool_BL7_10_14dpa   — 7-14 dpa  (GSM5045047, pooled)
#     Xen_BL14dpa             — 14 dpa    (GSM5057665, standalone)
#     Xen_Pool_BL14_20_52dpa  — 14-52 dpa (GSM5057660, pooled)
#
#   LBst (limb bud stage, REGENERATIVE tadpole):
#     Xen_LBst50              — NF stage 50 (GSM5057657)
#     Xen_LBst51              — NF stage 51 (GSM5057658)
#     Xen_LBst52              — NF stage 52 (GSM5057659)
#
# Metadata columns added here:
#   condition : "BL" (froglet, non-regenerative) | "LBst" (tadpole, regenerative)
#   timepoint : "0dpa","3dpa","7-14dpa","14dpa","14-52dpa" for BL;
#               "NF50","NF51","NF52" for LBst
#   sample    : individual library name (used by Harmony for batch correction)
#   species   : "Xenopus"
#
# Excluded (see 00_download_data.R for rationale):
#   GSM5057666  Xen_Pool_LBst54_BL0dpa  — mixed regen + non-regen library
#   GSM5057667  Xen_Transplant           — transplant experiment
# =============================================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# Ordered factor levels — BL timepoints first (dpa), then LBst (NF stage)
BL_ORDER   <- c("0dpa","3dpa","7-14dpa","14dpa","14-52dpa")
LBST_ORDER <- c("NF50","NF51","NF52")
TIMEPOINT_ORDER <- c(BL_ORDER, LBST_ORDER)

TIMEPOINT_COLORS <- c(
  "0dpa"      = "#37474F",  # BL: dark grey-blue
  "3dpa"      = "#5C6BC0",  # BL: indigo
  "7-14dpa"   = "#1565C0",  # BL: dark blue
  "14dpa"     = "#6A1B9A",  # BL: dark purple
  "14-52dpa"  = "#AD1457",  # BL: dark pink
  "NF50"      = "#1B5E20",  # LBst: dark green (regenerative)
  "NF51"      = "#388E3C",  # LBst: medium green
  "NF52"      = "#66BB6A"   # LBst: light green
)

# =============================================================================
# 1. LOAD RAW 10x DATA
# =============================================================================

load_10x <- function(path, project_name) {
  counts <- Read10X(data.dir = path)
  CreateSeuratObject(counts = counts, project = project_name, min.cells = 3)
}

raw <- file.path("data", "raw")
message("Loading Xenopus 10x data...")

# BL samples
xen_BL0      <- load_10x(file.path(raw, "Xen_BL0dpa"),             "Xen_BL0dpa")
xen_FACS_BL0 <- load_10x(file.path(raw, "Xen_FACS_BL0dpa"),        "Xen_FACS_BL0dpa")
xen_BL3      <- load_10x(file.path(raw, "Xen_BL3dpa"),             "Xen_BL3dpa")
xen_FACS_BL3 <- load_10x(file.path(raw, "Xen_FACS_BL3dpa"),        "Xen_FACS_BL3dpa")
xen_pool     <- load_10x(file.path(raw, "Xen_Pool_BL7_10_14dpa"),  "Xen_Pool_BL7_10_14dpa")
xen_14dpa    <- load_10x(file.path(raw, "Xen_BL14dpa"),            "Xen_BL14dpa")
xen_late     <- load_10x(file.path(raw, "Xen_Pool_BL14_20_52dpa"), "Xen_Pool_BL14_20_52dpa")

# LBst samples
xen_LBst50 <- load_10x(file.path(raw, "Xen_LBst50"), "Xen_LBst50")
xen_LBst51 <- load_10x(file.path(raw, "Xen_LBst51"), "Xen_LBst51")
xen_LBst52 <- load_10x(file.path(raw, "Xen_LBst52"), "Xen_LBst52")

all_objs <- list(xen_BL0, xen_FACS_BL0, xen_BL3, xen_FACS_BL3, xen_pool,
                 xen_14dpa, xen_late, xen_LBst50, xen_LBst51, xen_LBst52)
message(sprintf("Raw cells per sample: %s",
        paste(sapply(all_objs, ncol), collapse = " | ")))

# =============================================================================
# 2. PER-SAMPLE QC
# =============================================================================

add_mito <- function(seu) {
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

all_objs <- lapply(all_objs, add_mito)
list2env(setNames(all_objs,
  c("xen_BL0","xen_FACS_BL0","xen_BL3","xen_FACS_BL3","xen_pool",
    "xen_14dpa","xen_late","xen_LBst50","xen_LBst51","xen_LBst52")),
  envir = .GlobalEnv)

# QC violin — inspect before adjusting thresholds
p_qc <- VlnPlot(
  merge(xen_BL0, list(xen_FACS_BL0, xen_BL3, xen_FACS_BL3, xen_pool,
                      xen_14dpa, xen_late, xen_LBst50, xen_LBst51, xen_LBst52)),
  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
  ncol = 3, pt.size = 0
)
pdf(file.path(results_dir, "xenopus_qc_violin_prefilter.pdf"), width = 24, height = 4)
print(p_qc)
dev.off()
message("QC violin saved — inspect xenopus_qc_violin_prefilter.pdf before adjusting thresholds")

# BL thresholds: conservative starting point — inspect violin and adjust
xen_BL0      <- subset(xen_BL0,      nFeature_RNA > 500 & nFeature_RNA < 8000  & percent.mt < 20)
xen_FACS_BL0 <- subset(xen_FACS_BL0, nFeature_RNA > 500 & nFeature_RNA < 8000  & percent.mt < 20)
xen_BL3      <- subset(xen_BL3,      nFeature_RNA > 500 & nFeature_RNA < 12000 & percent.mt < 20)
xen_FACS_BL3 <- subset(xen_FACS_BL3, nFeature_RNA > 500 & nFeature_RNA < 12000 & percent.mt < 20)
xen_pool     <- subset(xen_pool,     nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)
xen_14dpa    <- subset(xen_14dpa,    nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)
xen_late     <- subset(xen_late,     nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)

# LBst thresholds: same conservative window — inspect violin and adjust
xen_LBst50 <- subset(xen_LBst50, nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)
xen_LBst51 <- subset(xen_LBst51, nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)
xen_LBst52 <- subset(xen_LBst52, nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)

message("Post-QC cells:")
message(sprintf(
  "  BL  : BL0=%d | FACS_BL0=%d | BL3=%d | FACS_BL3=%d | pool=%d | 14dpa=%d | late=%d",
  ncol(xen_BL0), ncol(xen_FACS_BL0), ncol(xen_BL3), ncol(xen_FACS_BL3),
  ncol(xen_pool), ncol(xen_14dpa), ncol(xen_late)
))
message(sprintf("  LBst: NF50=%d | NF51=%d | NF52=%d",
  ncol(xen_LBst50), ncol(xen_LBst51), ncol(xen_LBst52)
))

# =============================================================================
# 3. ADD METADATA
# =============================================================================

# condition: "BL" = froglet blastema (non-regenerative)
#            "LBst" = tadpole limb bud stage (regenerative)
# timepoint: dpa label for BL; NF stage for LBst

set_meta <- function(seu, condition, timepoint, sample) {
  seu$condition <- condition
  seu$timepoint <- timepoint
  seu$sample    <- sample
  seu$species   <- "Xenopus"
  seu
}

xen_BL0      <- set_meta(xen_BL0,      "BL",   "0dpa",     "Xen_BL0dpa")
xen_FACS_BL0 <- set_meta(xen_FACS_BL0, "BL",   "0dpa",     "Xen_FACS_BL0dpa")
xen_BL3      <- set_meta(xen_BL3,      "BL",   "3dpa",     "Xen_BL3dpa")
xen_FACS_BL3 <- set_meta(xen_FACS_BL3, "BL",   "3dpa",     "Xen_FACS_BL3dpa")
xen_pool     <- set_meta(xen_pool,     "BL",   "7-14dpa",  "Xen_Pool_BL7_10_14dpa")
xen_14dpa    <- set_meta(xen_14dpa,    "BL",   "14dpa",    "Xen_BL14dpa")
xen_late     <- set_meta(xen_late,     "BL",   "14-52dpa", "Xen_Pool_BL14_20_52dpa")
xen_LBst50   <- set_meta(xen_LBst50,  "LBst", "NF50",     "Xen_LBst50")
xen_LBst51   <- set_meta(xen_LBst51,  "LBst", "NF51",     "Xen_LBst51")
xen_LBst52   <- set_meta(xen_LBst52,  "LBst", "NF52",     "Xen_LBst52")

# =============================================================================
# 4. MERGE + HARMONY INTEGRATION
#
# 10 samples across two GEO submissions and two experimental conditions.
# Harmony corrects for library-level batch effects (group.by.vars = "sample").
# The condition and timepoint metadata survive batch correction — Harmony
# only adjusts the PCA embedding, not the expression values.
# =============================================================================

message("Merging all Xenopus objects...")
xen_merged <- merge(
  xen_BL0,
  y = list(xen_FACS_BL0, xen_BL3, xen_FACS_BL3, xen_pool,
           xen_14dpa, xen_late, xen_LBst50, xen_LBst51, xen_LBst52),
  add.cell.ids = c("BL0","FACS_BL0","BL3","FACS_BL3","pool",
                   "14dpa","late","LBst50","LBst51","LBst52"),
  project = "XenopusRegen"
)

# Free individual objects immediately to reclaim RAM before heavy steps
rm(xen_BL0, xen_FACS_BL0, xen_BL3, xen_FACS_BL3, xen_pool,
   xen_14dpa, xen_late, xen_LBst50, xen_LBst51, xen_LBst52, all_objs)
gc()

xen_merged$timepoint <- factor(xen_merged$timepoint, levels = TIMEPOINT_ORDER)
xen_merged$condition <- factor(xen_merged$condition, levels = c("BL","LBst"))

message("Total merged cells: ", ncol(xen_merged))
print(table(xen_merged$condition, xen_merged$timepoint))

# =============================================================================
# 5. NORMALISE, HVG, SCALE, PCA
# =============================================================================

message("Normalising and running PCA...")
xen_merged <- NormalizeData(xen_merged)
xen_merged <- FindVariableFeatures(xen_merged, nfeatures = 2000)
xen_merged <- ScaleData(xen_merged, features = VariableFeatures(xen_merged))
xen_merged <- RunPCA(xen_merged, npcs = 30)

pdf(file.path(results_dir, "xenopus_elbow_plot.pdf"), width = 6, height = 4)
print(ElbowPlot(xen_merged, ndims = 30))
dev.off()

gc()  # free dense scale.data matrix before Harmony

message("Running Harmony (batch correction by sample)...")
xen_merged <- RunHarmony(
  xen_merged,
  group.by.vars  = "sample",
  reduction      = "pca",
  reduction.save = "harmony",
  verbose        = FALSE
)
# JoinLayers is intentionally deferred to 02_xenopus_annotation.R so that
# script 01 never needs to hold split + joined layers simultaneously.

# =============================================================================
# 6. CLUSTERING + UMAP
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
print(table(xen_merged$condition, xen_merged$timepoint))

# UMAP by cluster
pdf(file.path(results_dir, "xenopus_umap_clusters.pdf"), width = 8, height = 6)
print(DimPlot(xen_merged, reduction = "umap_harmony", label = TRUE, repel = TRUE) +
        ggtitle("Xenopus - Seurat clusters (res=0.3)"))
dev.off()

# UMAP by timepoint (BL shades of blue/purple; LBst shades of green)
pdf(file.path(results_dir, "xenopus_umap_timepoint.pdf"), width = 10, height = 6)
print(DimPlot(xen_merged, reduction = "umap_harmony", group.by = "timepoint",
              cols = TIMEPOINT_COLORS, pt.size = 0.4) +
        ggtitle("Xenopus - by timepoint"))
dev.off()

# UMAP by condition — shows BL vs LBst separation
pdf(file.path(results_dir, "xenopus_umap_condition.pdf"), width = 8, height = 6)
print(DimPlot(xen_merged, reduction = "umap_harmony", group.by = "condition",
              cols = c("BL" = "#1565C0", "LBst" = "#2E7D32"), pt.size = 0.4) +
        ggtitle("Xenopus - BL (froglet, non-regen) vs LBst (tadpole, regen)"))
dev.off()

# Drop heavy layers before saving — counts and scale.data not needed after PCA/clustering
xen_merged[["RNA"]]$scale.data <- NULL
xen_merged[["RNA"]]$counts     <- NULL
gc()

saveRDS(xen_merged, file.path(results_dir, "xen_merged.rds"))
message("Saved results/xen_merged.rds")
