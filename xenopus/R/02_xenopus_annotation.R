# =============================================================================
# 02 — Xenopus laevis: cell type annotation
# Input:  results/xen_merged.rds
# Output: results/xen_merged.rds (updated with cell_type metadata)
#
# 20 clusters at resolution=0.3 (five timepoints: 0/3/7-14/14/14-52 dpa,
# Harmony-integrated across two GEO batches: GSM5045xxx + GSM5057xxx)
#
# Annotation strategy: marker-based verification via canonical dotplot.
# Each cluster is assigned a label only after inspecting its top expressed
# markers in the canonical_markers dotplot (Section 2). Labels are NOT
# assumed from cluster order. Where the dotplot evidence contradicts an
# initial assignment, the label is corrected to match marker identity.
#
# Key markers used:
#   col1a1.L/S, vim.L          — CT/stromal fibroblast
#   col2a1.L, sox9.L           — cartilage/chondrocyte
#   krt8.L, krt18.L, epcam.L  — epithelial/keratinocyte
#   hba1.L, hbb.S              — erythrocyte
#   ptprc.L                    — pan-immune (leukocyte)
#   pecam1.L, cdh5.L           — endothelial
#   acta2.L, tagln.L           — smooth muscle / myofibroblast
#   mki67.L, top2a.L           — proliferating
#   sox10.L                    — neural crest / Schwann
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
  "3"  = "Epithelial_mixed",      # epcam dominant — dotplot shows epithelial not CT (corrected)
  "4"  = "Smooth_muscle",         # col1a1, vim, acta2, tagln
  "12" = "Schwann_cell",          # sox10 dominant — dotplot shows neural/Schwann not tendon (corrected)
  "16" = "CT_cartilage_assoc",    # col2a1 — verify with hba1 panel (see canonical_markers below)
  # Epithelial
  "1"  = "Keratinocyte",          # epcam — strongest epithelial cluster
  "5"  = "Basal_epithelial",      # epcam, krt17, vim mix
  "8"  = "CT_progenitor",         # epcam, krt17 — progenitor-like mixed profile
  "10" = "Keratinocyte_2",        # krt17, epcam — top krt17 expression
  "18" = "Neural_2",              # sox10 dominant — dotplot shows neural not epithelial (corrected)
  "19" = "Keratinocyte_3",        # krt8, krt18 dominant — dotplot shows epithelial not osteoblast (corrected)
  # Immune
  "6"  = "Macrophage",            # ptprc — largest immune cluster
  "11" = "Neutrophil",            # ptprc, mki67 — cycling immune
  "17" = "Macrophage_2",          # ptprc — second macrophage subtype
  # Proliferating
  "7"  = "Proliferating_A",       # top2a, tagln — proliferating CT
  "9"  = "Proliferating_B",       # epcam, top2a — proliferating epithelial
  # Vascular / other
  "13" = "Endothelial",           # pecam1, cdh5 — only cluster with pecam1
  "14" = "Erythrocyte",           # hba1 — verify with hba1 panel (col1a1 signal needs check)
  "21" = "Erythrocyte_2",         # hba1 + krt18 — small late-stage cluster
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
  # Erythrocyte — included to verify clusters 14/16 (col1a1 signal in Erythrocyte row)
  "hba1.L","hbb.S",
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
  "Keratinocyte","Keratinocyte_2","Keratinocyte_3","Basal_epithelial","Epithelial_mixed",
  "CT_progenitor",
  "Macrophage","Macrophage_2","Neutrophil",
  "Endothelial","Smooth_muscle",
  "CT_fibroblast_A","CT_fibroblast_B",
  "CT_cartilage_assoc",
  "Proliferating_A","Proliferating_B",
  "Neural_crest","Neural_2","Schwann_cell","Neuron",
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
