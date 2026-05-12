# =============================================================================
# plot_hh_celltype_violin.R
#
# Violin plots of PTCH1, PTCH2, and GLI1 expression broken down by cell type
# across three species/datasets, for cross-species HH pathway comparison.
#
# Datasets:
#   Collaborator axolotl — annotated_celltype (from ShinyCell metadata)
#   Xenopus (Lin 2021)   — cell_type
#   Mouse digit tip      — seurat_clusters (not yet annotated)
#
# Run from comparison/ directory: source("R/plot_hh_celltype_violin.R")
#
# Output: results/hh_celltype_violin.pdf
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(hdf5r)

AXO_COLLAB_DIR <- "data/axo_collab"
XEN_RDS        <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
MOUSE_RDS      <- "~/Desktop/scrnaseq-pipelines/mouse-digit/results/mouse_merged.rds"
RESULTS_DIR    <- "results"

HH_GENES <- c("PTCH1", "PTCH2", "GLI1")

# =============================================================================
# HELPERS
# =============================================================================

violin_theme <- theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey93", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, color = "grey40"),
    legend.position    = "none"
  )

make_celltype_violin <- function(df, expr_col, celltype_col, gene_label,
                                 title, subtitle, ct_colors = NULL) {
  df2 <- df[!is.na(df[[celltype_col]]) & !is.na(df[[expr_col]]), ]
  # Order cell types by median expression descending
  med_order <- names(sort(tapply(df2[[expr_col]], df2[[celltype_col]], median),
                          decreasing = TRUE))
  df2[[celltype_col]] <- factor(df2[[celltype_col]], levels = med_order)

  n_ct <- length(unique(df2[[celltype_col]]))
  if (is.null(ct_colors)) {
    ct_colors <- setNames(
      colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                         "#FF7F00","#A65628","#F781BF","#999999"))(n_ct),
      levels(df2[[celltype_col]])
    )
  }

  ggplot(df2, aes(x = .data[[celltype_col]], y = .data[[expr_col]],
                  fill = .data[[celltype_col]])) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, color = "grey30") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",
                 color = "grey20", linewidth = 0.35) +
    scale_fill_manual(values = ct_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = title, subtitle = subtitle,
         x = NULL, y = paste0(gene_label, " (log-norm)")) +
    violin_theme
}

# Efficient extraction of a few genes from Seurat v5 split layers
fetch_genes_split <- function(seu, genes, meta_cols = character(0)) {
  all_layers  <- Layers(seu)
  data_layers <- all_layers[startsWith(all_layers, "data.") | all_layers == "data"]
  expr_list <- lapply(data_layers, function(ln) {
    mat  <- LayerData(seu, layer = ln)
    hits <- intersect(genes, rownames(mat))
    if (length(hits) == 0) return(NULL)
    df <- as.data.frame(t(as.matrix(mat[hits, , drop = FALSE])))
    for (g in setdiff(genes, hits)) df[[g]] <- 0
    df[, genes, drop = FALSE]
  })
  expr_df <- do.call(rbind, Filter(Negate(is.null), expr_list))
  if (length(meta_cols) > 0) {
    meta    <- seu@meta.data[rownames(expr_df), meta_cols, drop = FALSE]
    expr_df <- cbind(expr_df, meta)
  }
  expr_df
}

# =============================================================================
# 1. COLLABORATOR AXOLOTL — read from ShinyCell h5
# =============================================================================

message("Loading collaborator axolotl...")
axo_meta <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objmeta.rds"))
gene_idx  <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objgene.rds"))[["RNA"]]

# Get h5 row indices for our target genes
target_rows <- gene_idx[HH_GENES]
target_rows <- target_rows[!is.na(target_rows)]
found_genes <- names(target_rows)
message("  Found in axolotl: ", paste(found_genes, collapse=", "))

h5   <- H5File$new(file.path(AXO_COLLAB_DIR, "merged_objassay_RNA.h5"), mode = "r")
dset <- h5[["grp/data"]]

axo_expr <- setNames(
  lapply(found_genes, function(g) dset[target_rows[g], ]),
  found_genes
)
h5$close_all()

axo_df <- as.data.frame(axo_expr)
axo_df$cell_type <- axo_meta$annotated_celltype

n_axo <- nrow(axo_df)
message("  ", n_axo, " cells, ", length(unique(axo_df$cell_type)), " cell types")

# =============================================================================
# 2. XENOPUS — sum .L and .S homeologs
# =============================================================================

message("Loading Xenopus...")
xen <- readRDS(path.expand(XEN_RDS))

xen_gene_map <- list(
  PTCH1 = c("ptch1.L", "ptch1.S"),
  PTCH2 = c("ptch2.L", "ptch2.S"),
  GLI1  = c("gli1.L",  "gli1.S")
)
xen_fetch_genes <- unique(unlist(xen_gene_map))

xen_raw <- fetch_genes_split(xen, xen_fetch_genes, meta_cols = "cell_type")
rm(xen); gc()

xen_df <- data.frame(cell_type = xen_raw$cell_type)
for (gene in HH_GENES) {
  homeologs      <- intersect(xen_gene_map[[gene]], colnames(xen_raw))
  xen_df[[gene]] <- if (length(homeologs) > 0)
    rowSums(xen_raw[, homeologs, drop = FALSE]) else 0
}

message("  ", nrow(xen_df), " cells, ", length(unique(xen_df$cell_type)), " cell types")

# =============================================================================
# 3. MOUSE — use seurat_clusters (not yet annotated)
# =============================================================================

message("Loading Mouse...")
mus <- readRDS(path.expand(MOUSE_RDS))

mus_gene_map <- c(PTCH1 = "Ptch1", PTCH2 = "Ptch2", GLI1 = "Gli1")
mus_fetch    <- as.character(mus_gene_map)

mus_raw <- fetch_genes_split(mus, mus_fetch, meta_cols = "seurat_clusters")
rm(mus); gc()

mus_df <- data.frame(cell_type = paste0("C", mus_raw$seurat_clusters))
for (gene in HH_GENES) {
  mouse_sym      <- mus_gene_map[gene]
  mus_df[[gene]] <- if (mouse_sym %in% colnames(mus_raw))
    mus_raw[[mouse_sym]] else 0
}

message("  ", nrow(mus_df), " cells, ", length(unique(mus_df$cell_type)), " clusters")

# =============================================================================
# 4. BUILD PLOTS — one page per gene, three species per page
# =============================================================================

make_species_page <- function(gene) {
  p_axo <- make_celltype_violin(
    axo_df, gene, "cell_type", gene,
    title    = "Axolotl (collaborator)",
    subtitle = paste0(nrow(axo_df), " cells | regen + non-regen")
  )
  p_xen <- make_celltype_violin(
    xen_df, gene, "cell_type", gene,
    title    = "Xenopus (Lin 2021)",
    subtitle = paste0(nrow(xen_df), " cells | froglet + tadpole")
  )
  p_mus <- make_celltype_violin(
    mus_df, gene, "cell_type", gene,
    title    = "Mouse digit tip (Johnson 2020)",
    subtitle = paste0(nrow(mus_df), " cells | clusters (not yet annotated)")
  )

  (p_axo / p_xen / p_mus) +
    plot_annotation(
      title    = paste0(gene, " expression by cell type — cross-species comparison"),
      subtitle = "Cell types ordered by median expression (descending). Log-normalized values.",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40")
      )
    )
}

out_file <- file.path(RESULTS_DIR, "hh_celltype_violin.pdf")
pdf(out_file, width = 14, height = 16)
for (gene in HH_GENES) {
  message("Plotting ", gene, "...")
  print(make_species_page(gene))
}
dev.off()
message("Saved ", out_file)
