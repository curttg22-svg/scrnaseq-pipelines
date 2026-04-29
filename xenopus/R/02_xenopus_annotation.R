# =============================================================================
# 02 — Xenopus laevis: cell type annotation
# Input:  results/xen_merged.rds
# Output: results/xen_merged.rds (updated with cell_type metadata)
#
# 21 clusters at resolution=0.3 (10 samples: BL 0/3/7-14/14/14-52 dpa + LBst NF50/51/52,
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
Idents(xen_merged) <- xen_merged$seurat_clusters   # reset to cluster numbers

# =============================================================================
# 1. ASSIGN CELL TYPE LABELS
# Clusters determined at resolution=0.3 (21 clusters, 0-20)
# Evidence: canonical marker dotplot (see below)
# All labels re-verified after full 10-sample re-integration (BL + LBst)
# =============================================================================

xen_merged <- RenameIdents(xen_merged,
  # CT / stromal fibroblasts (col1a1/vim dominant; no acta2/tagln)
  "0"  = "CT_fibroblast_A",       # col1a1.L/S, vim.L strong
  "1"  = "CT_fibroblast_B",       # col1a1.L/S, vim.L strong
  "3"  = "CT_fibroblast_C",       # col1a1.L/S, vim.L
  "7"  = "CT_fibroblast_D",       # col1a1.L/S, vim.L
  "18" = "CT_fibroblast_E",       # col1a1.L/S, vim.L very strong (likely LBst CT)
  # Cartilage-associated
  "6"  = "CT_cartilage_assoc",    # col2a1.L, sox9.L dominant
  "15" = "CT_perichondrial",      # col1a1 + col2a1/sox9 moderate
  # Smooth muscle
  "13" = "Smooth_muscle",         # acta2.L, tagln.L very dominant
  # Epithelial
  "2"  = "Keratinocyte",          # epcam.L very dominant
  "5"  = "CT_progenitor",         # epcam moderate + col2a1/sox9 — progenitor-like
  "8"  = "Epithelial_mixed",      # krt8.L, epcam.L
  # Immune
  "10" = "Macrophage",            # ptprc.L very dominant
  "11" = "Macrophage_2",          # ptprc.L moderate
  # Vascular
  "12" = "Endothelial",           # pecam1.L, cdh5.L dominant
  # Erythrocyte
  "9"  = "Erythrocyte",           # hba1.L very dominant
  "14" = "Erythrocyte_2",         # hba1.L large + epcam moderate — re-verify
  "20" = "Erythrocyte_LBst",      # hba1.L dominant + col1a1/ptprc — flag for re-verify
  # Proliferating
  "4"  = "Proliferating_A",       # mki67.L, top2a.L + col1a1 — proliferating CT
  "17" = "Proliferating_B",       # mki67.L, top2a.L dominant
  # Neural / Schwann
  "16" = "Neural_crest",          # sox10.L very dominant
  "19" = "Schwann_cell"           # vim.L + sox10.L
)

xen_merged@meta.data[["cell_type"]] <- as.character(Idents(xen_merged))

# =============================================================================
# 2. CANONICAL MARKER DOTPLOT — evidence for annotations
# =============================================================================

canonical_markers <- c(
  # CT / stromal — both homeologs present; combined signal strengthens distinction
  "col1a1.L","col1a1.S","vim.L","vim.S",
  # Cartilage — both homeologs present
  "col2a1.L","col2a1.S","sox9.L","sox9.S",
  # Epithelial — both homeologs present
  "krt8.L","krt8.S","krt18.L","krt18.S","epcam.L","epcam.S",
  # Immune
  "ptprc.L",                        # macrophage/pan-immune (.S absent from object)
  "s100a8.L",                       # neutrophil (aif1/cd68/mpx absent from object)
  "cd79a.L",                        # B cell
  # Vascular — both homeologs present for pecam1
  "pecam1.L","pecam1.S","cdh5.L",
  # Muscle — both homeologs present
  "acta2.L","acta2.S","tagln.L","tagln.S",
  # Erythrocyte — both homeologs; hba1.S included to resolve cluster 14/16 ambiguity
  "hba1.L","hba1.S",
  # Proliferating — both where present (top2a.S absent from object)
  "mki67.L","mki67.S","top2a.L",
  # Neural — both homeologs present
  "sox10.L","sox10.S","s100b.L"
)

# Keep only markers present in the object
canonical_markers <- intersect(canonical_markers, rownames(xen_merged))

pdf(file.path(results_dir, "xenopus_canonical_marker_dotplot.pdf"),
    width = max(12, length(canonical_markers) * 0.5 + 3), height = 9)
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
  "Keratinocyte","Epithelial_mixed","CT_progenitor",
  "CT_fibroblast_A","CT_fibroblast_B","CT_fibroblast_C","CT_fibroblast_D","CT_fibroblast_E",
  "CT_cartilage_assoc","CT_perichondrial",
  "Smooth_muscle",
  "Macrophage","Macrophage_2",
  "Endothelial",
  "Erythrocyte","Erythrocyte_2","Erythrocyte_LBst",
  "Proliferating_A","Proliferating_B",
  "Neural_crest","Schwann_cell"
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

# =============================================================================
# 4. GOI MODULE SCORE ON UMAP
# Same gene list as 06_goi_visualization.R. Homeolog-aware lookup finds
# .L / .S variants present in the object so no manual suffix matching needed.
# =============================================================================

GOI <- c(
  "GLI1","PTCH1","PTCH2","SHH","IHH","DHH",
  "GDF5","MSX1","SALL1","GREM1",
  "FREM1","CYTL1","CLDN5","IL11","FGFBP1",
  "C3","LDLR","ALOX15B","MSMO1","CYP51A1","FDPS"
)

find_xen_genes <- function(symbols, obj) {
  all_genes <- rownames(obj)
  found <- c()
  for (g in symbols) {
    hits <- grep(paste0("^", g, "(\\.L|\\.S|\\.L\\.[0-9]+|\\.S\\.[0-9]+)?$"),
                 all_genes, ignore.case = TRUE, value = TRUE)
    found <- c(found, hits)
  }
  unique(found)
}

xen_goi <- find_xen_genes(GOI, xen_merged)
message("GOI homeologs found in object: ", length(xen_goi))

if (length(xen_goi) > 0) {
  xen_merged <- AddModuleScore(xen_merged, features = list(xen_goi),
                               name = "GOI_score")

  pdf(file.path(results_dir, "GOI_module_score_UMAP.pdf"), width = 8, height = 6)
  print(
    FeaturePlot(xen_merged, features = "GOI_score1",
                reduction = "umap_harmony", pt.size = 0.3, order = TRUE) +
      scale_color_gradient2(low = "lightgrey", mid = "#FDAE61",
                            high = "#B2182B", midpoint = 0) +
      ggtitle("Xenopus - GOI module score (Hedgehog + patterning)") +
      theme(legend.title = element_text(size = 8))
  )
  dev.off()

  # Violin: module score by cell type — shows which cell types drive the signal
  pdf(file.path(results_dir, "GOI_module_score_violin.pdf"), width = 12, height = 5)
  print(
    VlnPlot(xen_merged, features = "GOI_score1", group.by = "cell_type",
            pt.size = 0) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "none") +
      ggtitle("GOI module score by cell type")
  )
  dev.off()

  message("GOI module score figures saved")
} else {
  message("No GOI homeologs found in object — module score skipped")
}
