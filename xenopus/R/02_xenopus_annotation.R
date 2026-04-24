# =============================================================================
# 02 — Xenopus laevis: cell type annotation
# Input:  results/xen_merged.rds
# Output: results/xen_merged.rds (updated with cell_type metadata)
#
# 22 clusters at resolution=0.3 (five timepoints: 0/3/7-14/14/14-52 dpa)
# Cluster assignments based on canonical marker dotplot evidence (see below).
# Key markers used: epcam.L/krt17.L (epithelial), col2a1.L (cartilage),
#   ptprc.L/aif1.L (immune), pecam1.L (endothelial), hba1.L (erythrocyte),
#   sox10.L (neural), top2a.L/mki67.L (proliferating), acta2.L/tagln.L (muscle)
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
  # CT / stromal subtypes (col1a1/vim dominant; acta2/tagln distinguish muscle)
  "0"  = "CT_fibroblast_A",       # col1a1, vim, tagln, acta2
  "2"  = "CT_fibroblast_B",       # col1a1, vim, tagln, acta2
  "3"  = "CT_fibroblast_C",       # col1a1, vim — hepatocyte-like marker lost with new cells
  "4"  = "Smooth_muscle",         # col1a1, vim, acta2, tagln
  "12" = "Tendon_fibroblast",     # col1a1, vim, acta2, tagln
  "16" = "CT_cartilage_assoc",    # col2a1, col1a1 — only cluster with col2a1
  # Epithelial
  "1"  = "Keratinocyte",          # epcam, krt17 — strongest krt17 after cluster 10
  "5"  = "Basal_epithelial",      # epcam, krt17, vim mix
  "8"  = "CT_progenitor",         # epcam, krt17 — progenitor-like mixed profile
  "10" = "Keratinocyte_2",        # krt17, epcam — top krt17 expression
  "18" = "Specialized_epithelial",# epcam secondary
  "19" = "Osteoblast",            # small epcam cluster — retains osteoblast-like profile
  # Immune
  "6"  = "Macrophage",            # ptprc — largest immune cluster
  "11" = "Neutrophil",            # ptprc, mki67 — cycling immune
  "17" = "Macrophage_2",          # ptprc — second macrophage subtype
  # Proliferating
  "7"  = "Proliferating_A",       # top2a, tagln — proliferating CT
  "9"  = "Proliferating_B",       # epcam, top2a — proliferating epithelial
  # Vascular / other
  "13" = "Endothelial",           # pecam1 — only cluster with pecam1
  "14" = "Erythrocyte",           # hba1 — dominant erythrocyte cluster
  "21" = "Erythrocyte_2",         # hba1 + krt18 — small late-stage erythrocyte cluster
  # Neural
  "15" = "Neural_crest",          # sox10 — larger neural cluster
  "20" = "Neuron"                 # sox10 — smaller neural cluster (new timepoints)
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
    ggtitle("Xenopus - canonical markers by cell type")
)
dev.off()

# =============================================================================
# 3. ANNOTATED UMAP
# =============================================================================

CT_ORDER <- c(
  "Keratinocyte","Keratinocyte_2","Basal_epithelial","Specialized_epithelial",
  "CT_progenitor","Osteoblast",
  "Macrophage","Macrophage_2","Neutrophil",
  "Endothelial","Smooth_muscle",
  "CT_fibroblast_A","CT_fibroblast_B","CT_fibroblast_C","Tendon_fibroblast",
  "CT_cartilage_assoc",
  "Proliferating_A","Proliferating_B",
  "Neural_crest","Neuron",
  "Erythrocyte","Erythrocyte_2"
)
xen_merged$cell_type <- factor(as.character(xen_merged$cell_type), levels = CT_ORDER)
Idents(xen_merged)   <- xen_merged$cell_type

pdf(file.path(results_dir, "xenopus_umap_annotated.pdf"), width = 10, height = 7)
print(DimPlot(xen_merged, group.by = "cell_type", label = TRUE,
              repel = TRUE, pt.size = 0.3) +
        ggtitle("Xenopus laevis - limb regeneration blastema") +
        theme(legend.text = element_text(size = 8)))
dev.off()

saveRDS(xen_merged, file.path(results_dir, "xen_merged.rds"))
message("Saved results/xen_merged.rds with cell_type metadata")
message("Cell type counts:")
print(table(xen_merged$cell_type))
