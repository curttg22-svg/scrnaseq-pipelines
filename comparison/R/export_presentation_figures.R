# =============================================================================
# export_presentation_figures.R
#
# Exports all figures needed for the lab meeting presentation as PNG files
# in comparison/figures/. Renders axolotl collab UMAP plots directly and
# converts existing pipeline PDFs via pdftoppm.
#
# Run from comparison/ directory: source("R/export_presentation_figures.R")
# =============================================================================

library(ggplot2)
library(patchwork)
library(hdf5r)

AXO_COLLAB_DIR <- "data/axo_collab"
FIG_DIR        <- "figures"
XEN_RESULTS    <- "~/Desktop/scrnaseq-pipelines/xenopus/results"
MUS_RESULTS    <- "~/Desktop/scrnaseq-pipelines/mouse-digit/results"
COMP_RESULTS   <- "results"

dir.create(FIG_DIR, showWarnings = FALSE)

pdf_to_png <- function(pdf_path, out_png, page = 1, dpi = 200) {
  pdf_path <- path.expand(pdf_path)
  if (!file.exists(pdf_path)) { message("  SKIP (not found): ", pdf_path); return(invisible(NULL)) }
  tmp_prefix <- tempfile()
  cmd <- sprintf('pdftoppm -png -r %d -f %d -l %d "%s" "%s"',
                 dpi, page, page, pdf_path, tmp_prefix)
  system(cmd, ignore.stdout = TRUE)
  tmp_file <- sprintf("%s-%d.png", tmp_prefix, page)
  if (file.exists(tmp_file)) {
    file.rename(tmp_file, out_png)
    message("  Exported: ", basename(out_png))
  } else {
    message("  FAILED: ", basename(out_png))
  }
}

# =============================================================================
# 1. AXOLOTL COLLAB — render UMAP plots directly to PNG
# =============================================================================

message("Rendering axolotl collab UMAP plots...")

meta     <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objmeta.rds"))
dimr     <- readRDS(file.path(AXO_COLLAB_DIR, "merged_objdimr.rds"))
umap_mat <- dimr[["umap.cca"]]

umap_df <- data.frame(
  UMAP1     = umap_mat[, 1],
  UMAP2     = umap_mat[, 2],
  cell_type = meta$annotated_celltype,
  condition = factor(ifelse(grepl("^non-regenerating", meta$sample_id),
                            "NonRegen", "Regen"),
                     levels = c("Regen", "NonRegen"))
)

n_ct  <- length(unique(umap_df$cell_type))
ct_pal <- setNames(
  colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                     "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
                     "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"))(n_ct),
  sort(unique(umap_df$cell_type))
)

umap_base <- theme_void(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle    = element_text(size = 9, color = "grey40", hjust = 0.5),
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.key.height = unit(0.4, "cm"),
    plot.background  = element_rect(fill = "white", color = NA)
  )

p_ct <- ggplot(umap_df[order(is.na(umap_df$cell_type)), ],
               aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 0.2, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = ct_pal, name = "Cell type",
                     guide = guide_legend(override.aes = list(size = 3, alpha = 1),
                                          ncol = 1)) +
  labs(title = "Axolotl (collaborator) — cell types",
       subtitle = paste0(nrow(umap_df), " cells  |  ", n_ct, " cell types")) +
  umap_base

p_cond <- ggplot(umap_df[order(umap_df$condition), ],
                 aes(UMAP1, UMAP2, color = condition)) +
  geom_point(size = 0.2, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = c(Regen = "#1565C0", NonRegen = "#B71C1C"),
                     name = "Condition",
                     guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "Axolotl (collaborator) — condition",
       subtitle = "Regen (dpa3/14/23) vs NonRegen (limb)") +
  umap_base

ggsave(file.path(FIG_DIR, "axo_collab_umap_celltypes.png"),
       p_ct, width = 7, height = 5.5, dpi = 200, bg = "white")
message("  Exported: axo_collab_umap_celltypes.png")

ggsave(file.path(FIG_DIR, "axo_collab_umap_condition.png"),
       p_cond, width = 5.5, height = 4.5, dpi = 200, bg = "white")
message("  Exported: axo_collab_umap_condition.png")

rm(meta, dimr, umap_mat, umap_df, p_ct, p_cond); gc()

# =============================================================================
# 2. AXOLOTL COLLAB — housekeeping QC violin (page 1 = actin page)
# =============================================================================

message("Exporting housekeeping violin (actin page)...")
pdf_to_png(file.path(COMP_RESULTS, "housekeeping_violin_v2.pdf"),
           file.path(FIG_DIR, "axo_collab_housekeeping_violin.png"),
           page = 1, dpi = 180)

# =============================================================================
# 3. XENOPUS — convert existing PDFs
# =============================================================================

message("Exporting Xenopus figures...")
pdf_to_png(file.path(XEN_RESULTS, "xenopus_umap_timepoint.pdf"),
           file.path(FIG_DIR, "xenopus_umap_timepoint.png"), dpi = 200)
pdf_to_png(file.path(XEN_RESULTS, "xenopus_umap_clusters.pdf"),
           file.path(FIG_DIR, "xenopus_umap_clusters.png"), dpi = 200)
pdf_to_png(file.path(XEN_RESULTS, "xenopus_umap_condition.pdf"),
           file.path(FIG_DIR, "xenopus_umap_condition.png"), dpi = 200)

# =============================================================================
# 4. MOUSE — convert existing PDFs
# =============================================================================

message("Exporting Mouse figures...")
pdf_to_png(file.path(MUS_RESULTS, "mouse_umap_timepoint.pdf"),
           file.path(FIG_DIR, "mouse_umap_timepoint.png"), dpi = 200)

# =============================================================================
# 5. HH CELL TYPE VIOLIN — page 1 (PTCH1) from condition split
# =============================================================================

message("Exporting HH violin (condition split, page 1)...")
pdf_to_png(file.path(COMP_RESULTS, "hh_celltype_violin_by_condition.pdf"),
           file.path(FIG_DIR, "axo_collab_hh_violin_condition.png"),
           page = 4, dpi = 180)   # page 4 = summary dot plot

# GOI UMAP grid (page 2 = GOI_ALL)
message("Exporting GOI UMAP grid...")
pdf_to_png(file.path(COMP_RESULTS, "goi_umap_axo_collab.pdf"),
           file.path(FIG_DIR, "axo_collab_goi_umap.png"),
           page = 2, dpi = 180)

message("\nDone. Figures in ", FIG_DIR, "/:")
cat(paste(sort(list.files(FIG_DIR)), collapse = "\n"), "\n")
