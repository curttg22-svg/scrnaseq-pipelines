# =============================================================================
# build_cross_species_atlas.R
#
# Integrates axolotl, Xenopus, and mouse digit-tip scRNA-seq into a shared
# atlas using Seurat v5 + Harmony.
#
# Memory strategy:
#   Each species is processed independently, saved to disk, and freed before
#   the next is loaded. The merge step only loads the (much smaller) processed
#   objects. Set SUBSAMPLE_N = NULL to use all cells.
#
# Run from comparison/ directory: source("R/build_cross_species_atlas.R")
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

if (!requireNamespace("harmony", quietly = TRUE))
  stop("harmony required: install.packages('harmony')")
if (!requireNamespace("hdf5r", quietly = TRUE))
  stop("hdf5r required: install.packages('hdf5r')")

AXO_DIR    <- "data/axo_collab"
XEN_RDS    <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
MUS_RDS    <- "~/Desktop/scrnaseq-pipelines/mouse-digit/results/mouse_merged.rds"
ATLAS_RDS  <- "results/cross_species_atlas.rds"
FIG_DIR    <- "figures/atlas"
STAGE_DIR  <- "results/atlas_stage"   # temp storage for per-species objects

SUBSAMPLE_N    <- 20000   # cells per species (NULL = all cells)
N_HVG          <- 3000
N_PCS          <- 30
HARMONY_LAMBDA <- 1
CLUSTER_RES    <- 0.5
set.seed(42)

dir.create(FIG_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(STAGE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# =============================================================================
# 1. Determine 3-way common gene set
# =============================================================================
message("=== Step 1: Common gene set ===")

axo_gene_idx <- readRDS(file.path(AXO_DIR, "merged_objgene.rds"))[["RNA"]]
axo_hgnc     <- names(axo_gene_idx)

message("  Loading Xenopus gene names...")
xen_tmp  <- readRDS(path.expand(XEN_RDS))
xen_hgnc <- unique(toupper(sub("\\.(L|S)(\\.[0-9]+)?$", "", rownames(xen_tmp))))
rm(xen_tmp); gc(verbose = FALSE)

message("  Loading Mouse gene names...")
mus_tmp  <- readRDS(path.expand(MUS_RDS))
mus_hgnc <- unique(toupper(rownames(mus_tmp)))
rm(mus_tmp); gc(verbose = FALSE)

atlas_genes <- Reduce(intersect, list(axo_hgnc, xen_hgnc, mus_hgnc))
message(sprintf("  Axolotl %d | Xenopus %d | Mouse %d => 3-way: %d genes",
                length(axo_hgnc), length(xen_hgnc), length(mus_hgnc),
                length(atlas_genes)))

# =============================================================================
# 2. Axolotl: read from HDF5, subsample, build Seurat object, save to disk
# =============================================================================
axo_stage <- file.path(STAGE_DIR, "axo.rds")

if (!file.exists(axo_stage)) {
  message("=== Step 2: Axolotl ===")

  axo_meta <- as.data.frame(readRDS(file.path(AXO_DIR, "merged_objmeta.rds")))
  rownames(axo_meta) <- axo_meta$cellID

  # Subsample cell indices
  n_axo <- nrow(axo_meta)
  cell_idx <- if (!is.null(SUBSAMPLE_N) && SUBSAMPLE_N < n_axo) {
    sort(sample.int(n_axo, SUBSAMPLE_N))
  } else {
    seq_len(n_axo)
  }
  message(sprintf("  Using %d / %d cells", length(cell_idx), n_axo))

  # Gene order in h5: names(sort(gene_idx)) gives genes in h5 row order
  all_axo_genes  <- names(sort(axo_gene_idx))
  is_atlas_gene  <- all_axo_genes %in% atlas_genes
  atlas_in_order <- all_axo_genes[is_atlas_gene]
  n_axo_atlas    <- length(atlas_in_order)
  n_cells_use    <- length(cell_idx)

  message(sprintf("  Reading %d genes x %d cells from HDF5...",
                  n_axo_atlas, n_cells_use))

  # Accumulate as sparse triplets — avoids full dense matrix
  h5    <- hdf5r::H5File$new(file.path(AXO_DIR, "merged_objassay_RNA.h5"), mode = "r")
  dset  <- h5[["grp/data"]]
  chunk_size <- 1000
  n_all  <- length(all_axo_genes)
  n_chunks <- ceiling(n_all / chunk_size)

  tri_i <- list(); tri_j <- list(); tri_x <- list()
  dest_row <- 0L

  for (k in seq_len(n_chunks)) {
    g_start     <- (k - 1L) * chunk_size + 1L
    g_end       <- min(k * chunk_size, n_all)
    local_atlas <- which(is_atlas_gene[g_start:g_end])
    if (length(local_atlas) == 0L) next

    chunk <- dset[g_start:g_end, cell_idx]   # rows=chunk_genes, cols=sampled cells
    sub   <- chunk[local_atlas, , drop = FALSE]
    nz    <- which(sub != 0, arr.ind = TRUE)
    if (nrow(nz) > 0) {
      tri_i[[k]] <- nz[, 1L] + dest_row
      tri_j[[k]] <- nz[, 2L]
      tri_x[[k]] <- sub[nz]
    }
    dest_row <- dest_row + length(local_atlas)
    if (k %% 10 == 0) message("    chunk ", k, "/", n_chunks)
  }
  h5$close_all()

  axo_sparse <- sparseMatrix(
    i    = unlist(tri_i), j = unlist(tri_j), x = unlist(tri_x),
    dims = c(n_axo_atlas, n_cells_use),
    dimnames = list(atlas_in_order, axo_meta$cellID[cell_idx])
  )
  rm(tri_i, tri_j, tri_x); gc(verbose = FALSE)

  axo_meta_sub <- axo_meta[cell_idx, ]
  axo_meta_sub$species   <- "Axolotl"
  axo_meta_sub$condition <- ifelse(grepl("^non-regenerating", axo_meta_sub$sample_id),
                                   "NonRegen", "Regen")
  axo_meta_sub$cell_type <- as.character(axo_meta_sub$annotated_celltype)

  axo_seu <- CreateSeuratObject(
    counts    = axo_sparse,
    meta.data = axo_meta_sub[, c("species","condition","timepoint","cell_type","sample_id")],
    project   = "Axolotl"
  )
  axo_seu[["RNA"]]$data <- axo_sparse
  rm(axo_sparse); gc(verbose = FALSE)
  message(sprintf("  Axolotl: %d genes x %d cells", nrow(axo_seu), ncol(axo_seu)))

  saveRDS(axo_seu, axo_stage)
  message("  Saved: ", axo_stage)
  rm(axo_seu, axo_meta, axo_meta_sub); gc(verbose = FALSE)
} else {
  message("=== Step 2: Axolotl (cached) ===")
}

# =============================================================================
# 3. Xenopus: aggregate homeologs, subsample, build Seurat object, save
# =============================================================================
xen_stage <- file.path(STAGE_DIR, "xen.rds")

if (!file.exists(xen_stage)) {
  message("=== Step 3: Xenopus ===")

  xen_obj <- readRDS(path.expand(XEN_RDS))
  xen_obj <- JoinLayers(xen_obj)
  message("  Layers joined")

  # Subsample before extracting data (saves memory during extraction)
  n_xen <- ncol(xen_obj)
  xen_keep <- if (!is.null(SUBSAMPLE_N) && SUBSAMPLE_N < n_xen) {
    sample.int(n_xen, SUBSAMPLE_N)
  } else {
    seq_len(n_xen)
  }
  message(sprintf("  Using %d / %d cells", length(xen_keep), n_xen))
  xen_obj <- xen_obj[, xen_keep]

  xen_logmat <- GetAssayData(xen_obj, layer = "data")
  xen_genes  <- rownames(xen_logmat)
  xen_base   <- toupper(sub("\\.(L|S)(\\.[0-9]+)?$", "", xen_genes))

  keep_mask  <- xen_base %in% atlas_genes
  xen_logmat <- xen_logmat[keep_mask, ]
  xen_base   <- xen_base[keep_mask]
  message(sprintf("  %d homeolog rows -> %d HGNC base genes",
                  sum(keep_mask), length(unique(xen_base))))

  # Indicator matrix for homeolog averaging
  unique_bases <- atlas_genes[atlas_genes %in% unique(xen_base)]
  base_idx     <- match(xen_base, unique_bases)
  indicator    <- sparseMatrix(
    i = base_idx, j = seq_len(nrow(xen_logmat)), x = 1,
    dims = c(length(unique_bases), nrow(xen_logmat))
  )
  xen_agg  <- indicator %*% xen_logmat
  hcounts  <- Matrix::rowSums(indicator)
  xen_agg  <- xen_agg / hcounts
  rownames(xen_agg) <- unique_bases
  rm(indicator, xen_logmat); gc(verbose = FALSE)

  xen_meta <- xen_obj@meta.data
  if (!"species"   %in% colnames(xen_meta)) xen_meta$species   <- "Xenopus"
  if (!"cell_type" %in% colnames(xen_meta) &&
      "cell_type"  %in% colnames(xen_obj@meta.data)) {
    xen_meta$cell_type <- xen_obj$cell_type
  }
  keep_cols <- intersect(c("species","condition","timepoint","cell_type","sample"),
                          colnames(xen_meta))

  xen_seu <- CreateSeuratObject(
    counts    = xen_agg,
    meta.data = xen_meta[, keep_cols],
    project   = "Xenopus"
  )
  xen_seu[["RNA"]]$data <- xen_agg
  rm(xen_obj, xen_agg); gc(verbose = FALSE)
  message(sprintf("  Xenopus: %d genes x %d cells", nrow(xen_seu), ncol(xen_seu)))

  saveRDS(xen_seu, xen_stage)
  message("  Saved: ", xen_stage)
  rm(xen_seu, xen_meta); gc(verbose = FALSE)
} else {
  message("=== Step 3: Xenopus (cached) ===")
}

# =============================================================================
# 4. Mouse: uppercase genes, subsample, build Seurat object, save
# =============================================================================
mus_stage <- file.path(STAGE_DIR, "mus.rds")

if (!file.exists(mus_stage)) {
  message("=== Step 4: Mouse ===")

  mus_obj    <- readRDS(path.expand(MUS_RDS))
  mus_obj    <- JoinLayers(mus_obj)

  n_mus <- ncol(mus_obj)
  mus_keep <- if (!is.null(SUBSAMPLE_N) && SUBSAMPLE_N < n_mus) {
    sample.int(n_mus, SUBSAMPLE_N)
  } else {
    seq_len(n_mus)
  }
  message(sprintf("  Using %d / %d cells", length(mus_keep), n_mus))
  mus_obj <- mus_obj[, mus_keep]

  mus_logmat <- GetAssayData(mus_obj, layer = "data")
  mus_syms   <- toupper(rownames(mus_logmat))
  rownames(mus_logmat) <- mus_syms

  is_dup <- duplicated(mus_syms)
  if (any(is_dup)) {
    message("  Removing ", sum(is_dup), " duplicate HGNC symbols")
    mus_logmat <- mus_logmat[!is_dup, ]
  }
  mus_logmat <- mus_logmat[rownames(mus_logmat) %in% atlas_genes, ]

  mus_meta <- mus_obj@meta.data
  if (!"species" %in% colnames(mus_meta)) mus_meta$species <- "Mouse"
  keep_cols <- intersect(c("species","condition","timepoint","cell_type","sample"),
                          colnames(mus_meta))
  keep_cols <- union(keep_cols, "species")

  mus_seu <- CreateSeuratObject(
    counts    = mus_logmat,
    meta.data = mus_meta[, keep_cols],
    project   = "Mouse"
  )
  mus_seu[["RNA"]]$data <- mus_logmat
  rm(mus_obj, mus_logmat); gc(verbose = FALSE)
  message(sprintf("  Mouse: %d genes x %d cells", nrow(mus_seu), ncol(mus_seu)))

  saveRDS(mus_seu, mus_stage)
  message("  Saved: ", mus_stage)
  rm(mus_seu, mus_meta); gc(verbose = FALSE)
} else {
  message("=== Step 4: Mouse (cached) ===")
}

# =============================================================================
# 5. Merge and Harmony integration
# =============================================================================
message("=== Step 5: Merge and integrate ===")

axo_seu <- readRDS(axo_stage)
xen_seu <- readRDS(xen_stage)
mus_seu <- readRDS(mus_stage)

# Align gene sets — each species may have slightly different atlas gene coverage
common_final <- Reduce(intersect, list(rownames(axo_seu), rownames(xen_seu), rownames(mus_seu)))
message(sprintf("  Final gene set after per-species subsetting: %d genes", length(common_final)))
axo_seu <- axo_seu[common_final, ]
xen_seu <- xen_seu[common_final, ]
mus_seu <- mus_seu[common_final, ]

atlas <- merge(
  axo_seu,
  y            = list(xen_seu, mus_seu),
  add.cell.ids = c("axo", "xen", "mus"),
  project      = "CrossSpeciesAtlas"
)
rm(axo_seu, xen_seu, mus_seu); gc(verbose = FALSE)
atlas <- JoinLayers(atlas)
message(sprintf("  Merged: %d genes x %d cells", nrow(atlas), ncol(atlas)))

# Ensure species column is populated
if (any(is.na(atlas$species)))
  atlas$species[is.na(atlas$species)] <- atlas$orig.ident[is.na(atlas$species)]

message("  FindVariableFeatures...")
atlas <- FindVariableFeatures(atlas, nfeatures = N_HVG, verbose = FALSE)

message("  ScaleData...")
atlas <- ScaleData(atlas, verbose = FALSE)

message("  RunPCA...")
atlas <- RunPCA(atlas, npcs = N_PCS, verbose = FALSE)

message("  RunHarmony (by species)...")
atlas <- harmony::RunHarmony(
  atlas,
  group.by.vars  = "species",
  reduction      = "pca",
  reduction.save = "harmony",
  lambda         = HARMONY_LAMBDA,
  verbose        = FALSE
)

message("  FindNeighbors + FindClusters...")
atlas <- FindNeighbors(atlas, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
atlas <- FindClusters(atlas,  resolution = CLUSTER_RES, verbose = FALSE)

message("  RunUMAP...")
atlas <- RunUMAP(atlas, reduction = "harmony", dims = 1:N_PCS,
                 reduction.name = "umap.harmony", verbose = FALSE)

n_clust <- length(levels(atlas$seurat_clusters))
message(sprintf("  Clusters: %d  |  Cells: %d", n_clust, ncol(atlas)))

# =============================================================================
# 6. Save atlas
# =============================================================================
saveRDS(atlas, ATLAS_RDS)
message("Saved: ", ATLAS_RDS)

# =============================================================================
# 7. Core UMAP figures
# =============================================================================
message("=== Step 7: Figures ===")

umap_theme <- theme_void(base_size = 11) +
  theme(
    plot.background   = element_rect(fill = "#0d0d1a", color = NA),
    panel.background  = element_rect(fill = "#0d0d1a", color = NA),
    plot.title        = element_text(color = "grey97", face = "bold", size = 13,
                                     hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle     = element_text(color = "grey55", size = 8.5, hjust = 0.5),
    legend.text       = element_text(color = "grey85", size = 8),
    legend.title      = element_text(color = "grey70", size = 9, face = "bold"),
    legend.key        = element_rect(fill = NA, color = NA),
    legend.background = element_rect(fill = "#0d0d1a", color = NA),
    plot.margin       = margin(12, 12, 12, 12)
  )

species_pal <- c(Axolotl = "#1565C0", Xenopus = "#00897B", Mouse = "#43A047")
pt_size     <- pmin(0.6, 30000 / ncol(atlas))

umap_df <- as.data.frame(Embeddings(atlas, "umap.harmony"))
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$species   <- atlas$species
umap_df$cluster   <- atlas$seurat_clusters
umap_df$condition <- atlas$condition
umap_df$cell_type <- atlas$cell_type

shuffle <- sample(nrow(umap_df))

# --- 7a. By species ---
p_species <- ggplot(umap_df[shuffle, ], aes(UMAP1, UMAP2, color = species)) +
  geom_point(size = pt_size, alpha = 0.5, stroke = 0) +
  scale_color_manual(values = species_pal, name = "Species") +
  labs(title    = "Cross-species atlas — by species",
       subtitle = sprintf("%d cells  |  %d genes  |  Harmony (species)",
                          ncol(atlas), nrow(atlas))) +
  umap_theme +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(FIG_DIR, "atlas_umap_species.png"),
       p_species, width = 8, height = 6, dpi = 200, bg = "#0d0d1a")
message("  Saved atlas_umap_species.png")

# --- 7b. By cluster ---
clust_pal <- setNames(
  colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                     "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
                     "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
                     "#B3B3B3"))(n_clust),
  levels(atlas$seurat_clusters)
)
p_cluster <- ggplot(umap_df[shuffle, ], aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = pt_size, alpha = 0.5, stroke = 0) +
  scale_color_manual(values = clust_pal, name = "Cluster") +
  labs(title    = "Cross-species atlas — by cluster",
       subtitle = sprintf("%d clusters  |  res %.1f", n_clust, CLUSTER_RES)) +
  umap_theme +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                               ncol = ceiling(n_clust / 20)))
ggsave(file.path(FIG_DIR, "atlas_umap_clusters.png"),
       p_cluster, width = 8, height = 6, dpi = 200, bg = "#0d0d1a")
message("  Saved atlas_umap_clusters.png")

# --- 7c. Species facet ---
p_split <- ggplot(umap_df[shuffle, ], aes(UMAP1, UMAP2, color = species)) +
  geom_point(data = transform(umap_df[shuffle, ], species = NULL),
             color = "grey25", size = pt_size * 0.5, alpha = 0.15, stroke = 0) +
  geom_point(size = pt_size * 1.2, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = species_pal, guide = "none") +
  facet_wrap(~ species) +
  labs(title    = "Cross-species atlas — species split",
       subtitle = "Grey = all other cells") +
  umap_theme +
  theme(strip.text       = element_text(color = "grey90", size = 10, face = "bold"),
        strip.background = element_rect(fill = "#1a1a3a", color = NA))
ggsave(file.path(FIG_DIR, "atlas_umap_species_split.png"),
       p_split, width = 14, height = 5, dpi = 200, bg = "#0d0d1a")
message("  Saved atlas_umap_species_split.png")

# --- 7d. Axolotl cell type overlay ---
has_ct <- !is.na(umap_df$cell_type) & nzchar(as.character(umap_df$cell_type))
if (sum(has_ct) > 500) {
  ct_df   <- umap_df[has_ct, ]
  n_ct    <- length(unique(ct_df$cell_type))
  ct_pal  <- setNames(
    colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                       "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
                       "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"))(n_ct),
    unique(ct_df$cell_type)
  )
  p_ct <- ggplot(umap_df[shuffle, ], aes(UMAP1, UMAP2)) +
    geom_point(color = "grey18", size = pt_size * 0.5, alpha = 0.18, stroke = 0) +
    geom_point(data = ct_df[sample(nrow(ct_df)), ],
               aes(color = cell_type), size = pt_size * 1.4, alpha = 0.7, stroke = 0) +
    scale_color_manual(values = ct_pal, name = "Cell type\n(axolotl)") +
    labs(title    = "Cross-species atlas — axolotl cell type annotations",
         subtitle = "Annotated axolotl cells highlighted; all others in dark grey") +
    umap_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
  ggsave(file.path(FIG_DIR, "atlas_umap_celltype_axolotl.png"),
         p_ct, width = 10, height = 7, dpi = 200, bg = "#0d0d1a")
  message("  Saved atlas_umap_celltype_axolotl.png")
}

# --- 7e. Condition ---
cond_vals <- unique(na.omit(umap_df$condition))
if (length(cond_vals) >= 2) {
  cond_pal <- c(Regen = "#1565C0", NonRegen = "#B71C1C")
  valid    <- !is.na(umap_df$condition)
  p_cond   <- ggplot(umap_df[valid, ][sample(sum(valid)), ],
                     aes(UMAP1, UMAP2, color = condition)) +
    geom_point(size = pt_size, alpha = 0.5, stroke = 0) +
    scale_color_manual(values = cond_pal[intersect(names(cond_pal), cond_vals)],
                       name = "Condition") +
    labs(title    = "Cross-species atlas — regeneration condition",
         subtitle = "Regen = actively regenerating; NonRegen = intact / non-regen") +
    umap_theme +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  ggsave(file.path(FIG_DIR, "atlas_umap_condition.png"),
         p_cond, width = 8, height = 6, dpi = 200, bg = "#0d0d1a")
  message("  Saved atlas_umap_condition.png")
}

message(sprintf("\n=== Done: %d cells | %d genes | %d clusters ===",
                ncol(atlas), nrow(atlas), n_clust))
message("Atlas: ", ATLAS_RDS)
message("Figures: ", FIG_DIR)
