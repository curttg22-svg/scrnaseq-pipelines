# =============================================================================
# 02 — Xenopus laevis: cell type annotation
# Input:  results/xen_merged.rds
# Output: results/xen_merged.rds (updated with cell_type metadata)
# =============================================================================

library(Seurat)
library(ggplot2)

results_dir <- "results"
xen_merged  <- readRDS(file.path(results_dir, "xen_merged.rds"))

# =============================================================================
# 1. ASSIGN CELL TYPE LABELS
# Clusters determined at resolution=0.3 (20 clusters, 0-19)
# Evidence: canonical marker dotplot (see below)
# =============================================================================

xen_merged <- RenameIdents(xen_merged,
  "0"  = "CT_fibroblast_A",
  "1"  = "CT_cartilage_assoc",
  "2"  = "CT_fibroblast_B",
  "3"  = "Hepatocyte_like",
  "4"  = "Basal_epithelial",
  "5"  = "Keratinocyte",
  "6"  = "Macrophage",
  "7"  = "Proliferating_A",
  "8"  = "CT_progenitor",
  "9"  = "Neutrophil",
  "10" = "B_cell",
  "11" = "Smooth_muscle",
  "12" = "Specialized_epithelial",
  "13" = "Proliferating_B",
  "14" = "Endothelial",
  "15" = "Tendon_fibroblast",
  "16" = "Erythrocyte",
  "17" = "Neural_crest",
  "18" = "Neuron",
  "19" = "Osteoblast"
)

xen_merged$cell_type <- Idents(xen_merged)

# =============================================================================
# 2. CANONICAL MARKER DOTPLOT — evidence for annotations
# =============================================================================

canonical_markers <- c(
  # CT / stromal
  "col1a1.L","col1a1.S","vim.L",
  # Cartilage
  "col2a1.L","sox9.L",
  # Epithelial
  "krt8.L","krt18.L","epcam.L",
  # Immune
  "aif1.L","cd68.L","ptprc.L",    # macrophage
  "mpx.L","s100a8.L",              # neutrophil
  "cd79a.L",                       # B cell
  # Vascular
  "pecam1.L","cdh5.L",
  # Muscle
  "acta2.L","tagln.L",
  # Proliferating
  "mki67.L","top2a.L",
  # Neural
  "sox10.L","s100b.L"
)

# Keep only markers present in the object
canonical_markers <- intersect(canonical_markers, rownames(xen_merged))

pdf(file.path(results_dir, "xenopus_canonical_marker_dotplot.pdf"),
    width = max(10, length(canonical_markers) * 0.55 + 3), height = 8)
print(
  DotPlot(xen_merged, features = canonical_markers, group.by = "cell_type",
          dot.scale = 6, col.min = 0) +
    scale_color_gradient(low = "lightgrey", high = "#08306B") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8)) +
    ggtitle("Xenopus — canonical markers by cell type")
)
dev.off()

# =============================================================================
# 3. ANNOTATED UMAP
# =============================================================================

CT_ORDER <- c(
  "Keratinocyte","Basal_epithelial","Specialized_epithelial",
  "Macrophage","Neutrophil","B_cell",
  "Endothelial","Smooth_muscle",
  "CT_fibroblast_A","CT_fibroblast_B","Tendon_fibroblast",
  "CT_progenitor","CT_cartilage_assoc","Osteoblast",
  "Hepatocyte_like","Proliferating_A","Proliferating_B",
  "Neural_crest","Neuron","Erythrocyte"
)
xen_merged$cell_type <- factor(as.character(xen_merged$cell_type), levels = CT_ORDER)
Idents(xen_merged)   <- xen_merged$cell_type

pdf(file.path(results_dir, "xenopus_umap_annotated.pdf"), width = 10, height = 7)
print(DimPlot(xen_merged, group.by = "cell_type", label = TRUE,
              repel = TRUE, pt.size = 0.3) +
        ggtitle("Xenopus laevis — limb regeneration blastema") +
        theme(legend.text = element_text(size = 8)))
dev.off()

saveRDS(xen_merged, file.path(results_dir, "xen_merged.rds"))
message("Saved results/xen_merged.rds with cell_type metadata")
message("Cell type counts:")
print(table(xen_merged$cell_type))
