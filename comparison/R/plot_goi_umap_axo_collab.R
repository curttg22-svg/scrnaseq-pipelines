# =============================================================================
# plot_goi_umap_axo_collab.R
#
# UMAP feature plots for genes of interest in the collaborator axolotl dataset.
# Reads UMAP coordinates from merged_objdimr.rds and expression values
# directly from the ShinyCell HDF5 (one row per gene — fast).
#
# Two pages of GOI panels + one page for HH signaling components.
# Cells colored by log-norm expression; missing genes noted in subtitle.
#
# Run from comparison/ directory:
#   source("R/plot_goi_umap_axo_collab.R")
#
# Output: results/goi_umap_axo_collab.pdf
# =============================================================================

library(ggplot2)
library(patchwork)
library(hdf5r)

AXO_COLLAB_DIR <- "data/axo_collab"
RESULTS_DIR    <- "results"

GOI_ALL <- c(
  "C3", "LDLR", "MSMO1", "CYP51A1", "FGFBP1", "FDPS",
  "GDF5", "IL11",  "MSX1",  "CLDN5",  "CYTL1",  "FREM1",
  "GLI1", "GREM1", "PTCH1", "PTCH2"
)

GOI_HH <- c(
  "SHH", "IHH", "DHH",
  "PTCH1", "PTCH2",
  "SMO",
  "GLI1", "GLI2", "GLI3",
  "HHIP",
  "BOC", "GAS1", "CDON"
)

# =============================================================================
# 1. LOAD UMAP COORDS + METADATA + GENE INDEX
# =============================================================================

message("Loading UMAP coordinates...")
dimr      <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objdimr.rds"))
umap_mat  <- dimr[["umap.cca"]]            # 43083 x 2
umap_df   <- data.frame(
  UMAP1 = umap_mat[, 1],
  UMAP2 = umap_mat[, 2]
)

message("Loading metadata and gene index...")
meta     <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objmeta.rds"))
gene_idx <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objgene.rds"))[["RNA"]]

umap_df$cell_type <- meta$annotated_celltype
umap_df$condition <- ifelse(grepl("^non-regenerating", meta$sample_id),
                            "NonRegen", "Regen")

message("  ", nrow(umap_df), " cells  |  ",
        length(unique(umap_df$cell_type)), " cell types")

# =============================================================================
# 2. FETCH SINGLE GENE EXPRESSION FROM H5
# =============================================================================

h5   <- H5File$new(file.path(AXO_COLLAB_DIR, "merged_objassay_RNA.h5"), mode = "r")
dset <- h5[["grp/data"]]

fetch_gene <- function(symbol) {
  row <- gene_idx[symbol]
  if (is.na(row)) return(NULL)
  dset[row, ]
}

# Pre-fetch all unique genes needed
all_genes  <- unique(c(GOI_ALL, GOI_HH))
message("Fetching ", length(all_genes), " genes from h5...")
expr_cache <- list()
for (g in all_genes) {
  v <- fetch_gene(g)
  if (!is.null(v)) {
    expr_cache[[g]] <- v
  } else {
    message("  Not found: ", g)
  }
}
h5$close_all()

found_all <- intersect(GOI_ALL, names(expr_cache))
found_hh  <- intersect(GOI_HH,  names(expr_cache))
message("GOI_ALL found: ", length(found_all), "/", length(GOI_ALL))
message("GOI_HH  found: ", length(found_hh), "/", length(GOI_HH))

# =============================================================================
# 3. PLOT HELPERS
# =============================================================================

umap_theme <- theme_void(base_size = 9) +
  theme(
    plot.title    = element_text(face = "bold.italic", size = 9, hjust = 0.5),
    plot.subtitle = element_text(size = 6, color = "grey50", hjust = 0.5),
    legend.position = "right",
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    legend.title  = element_text(size = 7),
    legend.text   = element_text(size = 6)
  )

# Color scale: white -> gold -> red (intuitive for expression)
EXPR_COLORS <- c("#E8E8E8", "#FFF176", "#FF8F00", "#B71C1C")

make_feature_plot <- function(symbol, point_size = 0.15) {
  if (!symbol %in% names(expr_cache)) {
    # Return a blank placeholder labelled "not detected"
    return(
      ggplot(umap_df, aes(UMAP1, UMAP2)) +
        geom_point(size = point_size, color = "#E0E0E0", alpha = 0.4) +
        labs(title = symbol, subtitle = "not detected") +
        umap_theme
    )
  }
  expr  <- expr_cache[[symbol]]
  cap   <- quantile(expr, 0.99)   # cap at 99th pct to avoid outlier dominance
  expr_c <- pmin(expr, cap)

  df <- umap_df
  df$expr <- expr_c
  # Plot expressing cells on top
  df <- df[order(df$expr), ]

  hi <- max(cap, 0.01)   # guard against all-zero gene
  ggplot(df, aes(UMAP1, UMAP2, color = expr)) +
    geom_point(size = point_size, alpha = 0.7, stroke = 0) +
    scale_color_gradientn(
      colors = EXPR_COLORS,
      name   = "log-norm",
      limits = c(0, hi)
    ) +
    labs(title = symbol,
         subtitle = sprintf("max=%.2f  99pct=%.2f", max(expr), cap)) +
    umap_theme
}

# Build a patchwork page from a vector of gene symbols
make_gene_page <- function(genes, ncol = 4, page_title = NULL, page_subtitle = NULL) {
  plots <- lapply(genes, make_feature_plot)
  pg <- wrap_plots(plots, ncol = ncol)
  if (!is.null(page_title)) {
    pg <- pg + plot_annotation(
      title    = page_title,
      subtitle = page_subtitle,
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40")
      )
    )
  }
  pg
}

# =============================================================================
# 4. CELL TYPE + CONDITION REFERENCE UMAPS
# =============================================================================

ct_pal <- setNames(
  colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                     "#FF7F00","#A65628","#F781BF","#999999",
                     "#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
                     "#A6D854","#FFD92F","#E5C494"))(length(unique(umap_df$cell_type))),
  sort(unique(umap_df$cell_type))
)

p_ct <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 0.15, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = ct_pal, name = "Cell type",
                     guide = guide_legend(override.aes = list(size = 2.5, alpha = 1),
                                          ncol = 1)) +
  labs(title = "Cell type", subtitle = "Collaborator axolotl") +
  umap_theme +
  theme(legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.35, "cm"))

p_cond <- ggplot(umap_df, aes(UMAP1, UMAP2, color = condition)) +
  geom_point(size = 0.15, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = c(Regen = "#1565C0", NonRegen = "#B71C1C"),
                     name = "Condition",
                     guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = "Condition", subtitle = "Regen vs NonRegen") +
  umap_theme +
  theme(legend.position = "right")

# =============================================================================
# 5. SAVE PDF
# =============================================================================

missing_all <- setdiff(GOI_ALL, names(expr_cache))
missing_hh  <- setdiff(GOI_HH,  names(expr_cache))

out_file <- file.path(RESULTS_DIR, "goi_umap_axo_collab.pdf")
pdf(out_file, width = 16, height = 14)

# Reference page
message("Printing reference UMAPs...")
print(
  (p_ct | p_cond) +
    plot_annotation(
      title    = "Collaborator axolotl UMAP — reference",
      subtitle = paste0(nrow(umap_df), " cells  |  ",
                        length(unique(umap_df$cell_type)), " cell types  |  ",
                        "Regen + NonRegen conditions"),
      theme = theme(plot.title = element_text(face = "bold", size = 13),
                    plot.subtitle = element_text(size = 9, color = "grey40"))
    )
)

# GOI_ALL — two rows of 8
message("Printing GOI_ALL feature plots...")
print(
  make_gene_page(
    GOI_ALL, ncol = 4,
    page_title    = "Genes of interest — UMAP expression (axolotl collaborator)",
    page_subtitle = paste0(
      "Log-normalized expression, capped at 99th percentile per gene.  |  ",
      if (length(missing_all) > 0)
        paste0("Not detected: ", paste(missing_all, collapse = ", "))
      else "All genes detected."
    )
  )
)

# GOI_HH — 13 genes, 4 per row
message("Printing GOI_HH feature plots...")
print(
  make_gene_page(
    GOI_HH, ncol = 4,
    page_title    = "Hedgehog signaling components — UMAP expression (axolotl collaborator)",
    page_subtitle = paste0(
      "Log-normalized expression, capped at 99th percentile per gene.  |  ",
      if (length(missing_hh) > 0)
        paste0("Not detected: ", paste(missing_hh, collapse = ", "))
      else "All genes detected."
    )
  )
)

dev.off()
message("Saved ", out_file)
