# =============================================================================
# 06 — GOI visualization: expression across cell types and timepoints
# Input:  results/xen_merged.rds, results/axo_annotated.rds,
#         results/final_bridge.rds
# Output: PDF figures in results/
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

results_dir <- "results"

xen_merged    <- readRDS(file.path(results_dir, "xen_merged.rds"))
axo_annotated <- readRDS(file.path(results_dir, "axo_annotated.rds"))
final_bridge  <- readRDS(file.path(results_dir, "final_bridge.rds"))

# =============================================================================
# 1. GENES OF INTEREST
# =============================================================================

GOI <- c(
  "GLI1","PTCH1","PTCH2","SHH","IHH","DHH",       # Hedgehog pathway
  "GDF5","MSX1","SALL1","GREM1",                   # Patterning / BMP
  "FREM1","CYTL1","CLDN5","IL11","FGFBP1",         # CT / vascular / immune
  "C3","LDLR","ALOX15B","MSMO1","CYP51A1","FDPS"  # Immune / lipid metabolism
)

# =============================================================================
# 2. XENOPUS — find homeolog-aware gene names (.L / .S)
# =============================================================================

find_xen_genes <- function(gene_symbols, xen_object) {
  xen_all <- rownames(xen_object)
  found   <- c()
  for (g in gene_symbols) {
    matches <- grep(
      paste0("^", g, "(\\.L|\\.S|\\.L\\.[0-9]+|\\.S\\.[0-9]+)?$"),
      xen_all, ignore.case = TRUE, value = TRUE
    )
    found <- c(found, matches)
  }
  unique(found)
}

xen_goi  <- find_xen_genes(GOI, xen_merged)
message("Xenopus GOI found: ", length(xen_goi), "/", length(GOI))

# =============================================================================
# 3. AXOLOTL — map GOI to AMEX IDs
# =============================================================================

# Bridge column detection: handle both naming conventions
# final_bridge from GTF pipeline uses amex_symbol; NCBI ortholog bridge uses human_symbol
sym_col <- if ("amex_symbol" %in% names(final_bridge)) "amex_symbol" else "human_symbol"
id_col  <- if ("amex_symbol" %in% names(final_bridge)) "amex_id"     else "amex_id"

axo_goi_map <- final_bridge |>
  filter(toupper(.data[[sym_col]]) %in% GOI,
         .data[[id_col]] %in% rownames(axo_annotated)) |>
  mutate(gene = toupper(.data[[sym_col]])) |>
  select(amex_id = all_of(id_col), gene) |>
  distinct()

message("Axolotl GOI found: ", length(unique(axo_goi_map$gene)), "/", length(GOI))
if (length(unique(axo_goi_map$gene)) == 0)
  message("  Note: bridge IDs do not match axolotl object rownames — axolotl dotplot skipped")

# =============================================================================
# 4. XENOPUS DOTPLOT — GOI x cell type
# =============================================================================

pdf(file.path(results_dir, "xenopus_goi_dotplot.pdf"),
    width = max(10, length(xen_goi) * 0.5 + 3), height = 9)
print(
  DotPlot(xen_merged, features = xen_goi, group.by = "cell_type",
          dot.scale = 6, col.min = 0) +
    scale_color_gradient(low = "lightgrey", high = "#08306B") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
          axis.text.y = element_text(size = 8)) +
    ggtitle("Xenopus - GOI expression by cell type")
)
dev.off()

# =============================================================================
# 5. XENOPUS DOTPLOT — GOI x timepoint
# =============================================================================

pdf(file.path(results_dir, "xenopus_goi_dotplot_timepoint.pdf"),
    width = max(10, length(xen_goi) * 0.5 + 3), height = 5)
print(
  DotPlot(xen_merged, features = xen_goi, group.by = "timepoint",
          dot.scale = 7, col.min = 0) +
    scale_color_gradient(low = "lightgrey", high = "#08306B") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    ggtitle("Xenopus - GOI expression by timepoint")
)
dev.off()

# =============================================================================
# 6. AXOLOTL DOTPLOT — GOI x cell type
# Safe sparse extractor — avoids 5 GB dense conversion
# =============================================================================

get_data_sparse <- function(seu, assay = "RNA") {
  tryCatch(
    GetAssayData(seu, assay = assay, layer = "data"),
    error = function(e) seu[[assay]]$data
  )
}

if (nrow(axo_goi_map) > 0) {
  axo_plot    <- axo_annotated
  mat_orig    <- get_data_sparse(axo_plot)
  renamed_rows <- rownames(mat_orig)
  for (i in seq_len(nrow(axo_goi_map))) {
    idx <- which(renamed_rows == axo_goi_map$amex_id[i])
    if (length(idx) == 1) renamed_rows[idx] <- axo_goi_map$gene[i]
  }
  rownames(mat_orig) <- renamed_rows
  axo_plot[["RNA_named"]] <- CreateAssay5Object(counts = mat_orig)
  DefaultAssay(axo_plot)  <- "RNA_named"

  axo_goi_present <- intersect(unique(axo_goi_map$gene),
                                rownames(axo_plot[["RNA_named"]]))
  pdf(file.path(results_dir, "axolotl_goi_dotplot.pdf"),
      width = max(8, length(axo_goi_present) * 0.55 + 3), height = 5)
  print(
    DotPlot(axo_plot, features = axo_goi_present, group.by = "cell_type_v2",
            assay = "RNA_named", dot.scale = 7, col.min = 0) +
      scale_color_gradient(low = "lightgrey", high = "#B2182B") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      ggtitle("Axolotl - GOI expression by cell type (11 dpa)")
  )
  dev.off()
} else {
  message("Axolotl dotplot skipped — no GOI matched in bridge")
}

message("GOI figures saved to: ", results_dir)
