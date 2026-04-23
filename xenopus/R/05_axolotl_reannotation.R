# =============================================================================
# 05 — Axolotl: v2 cell type annotation + marker discovery
# Input:  results/axo_final.rds, results/final_bridge.rds
# Output: results/axo_annotated.rds, results/axo_markers.rds
#
# The 7 clusters at res=0.5 are annotated using a canonical marker dotplot.
# Clusters 0,1,2 = CT blastema (fate-map reporter+)
# Clusters 3,5   = CT subpopulations
# Cluster 4      = Non-CT_A (Macrophage by AIF1/CD68/PTPRC)
# Cluster 6      = Non-CT_B (reporter-negative CT)
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

results_dir <- "results"
axo_final    <- readRDS(file.path(results_dir, "axo_final.rds"))
final_bridge <- readRDS(file.path(results_dir, "final_bridge.rds"))

# =============================================================================
# 1. BUILD AMEX -> DISPLAY LABEL MAP
# =============================================================================

# Genes to look up by human symbol
lookup_genes <- c(
  "COL1A1","COL2A1","VIM","FN1",          # CT / stromal
  "AIF1","CD68","PTPRC","CSF1R",           # macrophage
  "PECAM1","CDH5",                          # endothelial
  "ACTA2","TAGLN",                          # smooth muscle
  "KRT8","KRT18","EPCAM",                  # epithelial
  "MKI67","TOP2A",                          # proliferating
  "SOX10","S100B"                           # neural/Schwann
)

amex_lookup <- final_bridge |>
  filter(amex_symbol %in% toupper(lookup_genes)) |>
  filter(amex_id %in% rownames(axo_final)) |>
  distinct(amex_symbol, .keep_all = TRUE)

message("Canonical markers found in axolotl object: ",
        nrow(amex_lookup), "/", length(lookup_genes))

# =============================================================================
# 2. CANONICAL MARKER DOTPLOT
# =============================================================================

marker_ids <- amex_lookup$amex_id
names(marker_ids) <- amex_lookup$amex_symbol

pdf(file.path(results_dir, "axolotl_canonical_marker_dotplot.pdf"),
    width = max(10, length(marker_ids) * 0.6 + 3), height = 5)
print(
  DotPlot(axo_final, features = unname(marker_ids), dot.scale = 7, col.min = 0) +
    scale_color_gradient(low = "lightgrey", high = "#08306B") +
    scale_x_discrete(labels = names(marker_ids)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    ggtitle("Axolotl — canonical markers by cluster")
)
dev.off()

# =============================================================================
# 3. ORIGINAL CLUSTER LABELS (v1)
# =============================================================================

axo_final <- RenameIdents(axo_final,
  "0" = "CT_blastema_A",
  "1" = "CT_blastema_B",
  "2" = "CT_transitional",
  "3" = "CT_differentiating",
  "4" = "Non_CT_A",
  "5" = "CT_specialized",
  "6" = "Non_CT_B"
)
axo_final$cell_type <- Idents(axo_final)

# =============================================================================
# 4. REVISED LABELS (v2) — evidence from canonical marker dotplot
#    Non_CT_A = Macrophage (AIF1+, CD68+, PTPRC+, CSF1R+)
# =============================================================================

new_annotation <- c(
  "CT_blastema_A"      = "Blastema_CT_A",       # ~2,558 cells — COL1A1+, reporter+
  "CT_blastema_B"      = "Blastema_CT_B",       # ~2,287 cells — COL1A1+, reporter+
  "CT_transitional"    = "CT_progenitor",        #   ~142 cells — transitional markers
  "CT_differentiating" = "CT_differentiating",   #   ~141 cells — COL1A1+, KRT18+ mixed
  "CT_specialized"     = "CT_specialized",       #    ~85 cells — distinct CT profile
  "Non_CT_A"           = "Macrophage",           #   ~121 cells — AIF1+, CD68+, PTPRC+
  "Non_CT_B"           = "CT_unreported"         #    ~57 cells — COL1A1+, reporter-
)

# Seurat v5: use AddMetaData, not $ assignment (avoids "No cell overlap" error)
ct_v2_vec <- unname(new_annotation[as.character(axo_final$cell_type)])
axo_annotated <- AddMetaData(axo_final,
                               metadata = ct_v2_vec,
                               col.name = "cell_type_v2")
Idents(axo_annotated) <- axo_annotated$cell_type_v2

# =============================================================================
# 5. ANNOTATED UMAP
# =============================================================================

pdf(file.path(results_dir, "axolotl_umap_annotated.pdf"), width = 8, height = 6)
print(DimPlot(axo_annotated, group.by = "cell_type_v2",
              label = TRUE, repel = TRUE, pt.size = 0.5) +
        ggtitle("Axolotl — A. mexicanum 11 dpa (CT fate-map enriched)"))
dev.off()

message("Cell type v2 counts:")
print(table(axo_annotated$cell_type_v2))

# =============================================================================
# 6. FIND ALL MARKERS
# =============================================================================

message("Running FindAllMarkers (this may take several minutes)...")
axo_markers <- FindAllMarkers(axo_annotated,
                                only.pos = TRUE,
                                min.pct  = 0.25,
                                logfc.threshold = 0.25)

saveRDS(axo_markers,    file.path(results_dir, "axo_markers.rds"))
saveRDS(axo_annotated,  file.path(results_dir, "axo_annotated.rds"))
message("Saved results/axo_annotated.rds and results/axo_markers.rds")
