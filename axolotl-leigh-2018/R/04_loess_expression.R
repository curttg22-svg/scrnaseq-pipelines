# =============================================================================
# 04 — Temporal expression curves: GOI trajectory across regeneration stages
# Input:  results/axo_integrated.rds
# Output: results/expression_curves/<GENE>_expr_curve.pdf (one per GOI)
#
# For each GOI, plots mean log-normalized expression (left panel) and percent
# of cells expressing (right panel) across the four regeneration stages
# (Intact -> WH -> EB -> MB), with one line per cell type.
#
# Only cell types with detectable expression (mean > 0.05 in at least one
# timepoint) are shown per gene, preventing clutter from non-expressing types.
# WH_N4 / WH_N5 / WH_N6 are pooled into "WH" before computing means.
#
# Rationale for design choice:
#   Timepoints are biologically ordered (real regeneration stages), so
#   connecting them with lines shows a meaningful temporal trajectory for
#   each cell type. The previous design placed categorical cell types on
#   the x-axis and connected them with lines — that implied a continuous
#   relationship between unrelated cell types that does not exist, and was
#   replaced here. Dotplots in script 03 already handle cell-type-by-cell-type
#   comparisons; this script is for temporal trajectory only.
#
# Note: this script uses single-cell data only. Bulk RNA-seq comparisons
# are not included (private data, not redistributable).
# =============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

results_dir <- "results"
out_dir     <- file.path(results_dir, "expression_curves")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TIMEPOINT_ORDER <- c("Intact", "WH", "EB", "MB")
WH_REPS         <- c("WH_N4", "WH_N5", "WH_N6")

# Biologically ordered cell types (paper labels — unlisted types appended)
CT_ORDER <- c(
  "Macrophage","Recruited_Macrophage","Phagocytosing_Neutrophil",
  "Neutrophil_#1","Neutrophil_#2","T_cell","Early_B",
  "Erythrocyte","Erythrocyte #1","Erythrocyte #2",
  "Intermediate_epidermis","Proliferating_Epidermis","BasalEpidermis",
  "Epidermal_Langerhans","Basal wound epidermis","Basal WE ",
  "Intermediate_WE_#1","Intermediate_WE_#2","Intermediate_WE_#3",
  "Intermediate WE #4","Intermediate WE #5",
  "SSC","SSCs","Fibroblast",
  "Fibroblast-like_blastema_#1","Fibroblast-like_blastema_#2",
  "Fibroblast-like blastema #3","Fibroblast-like blastema #4",
  "Fibroblast-like blastema #5","Fibroblast-like blastema #7",
  "Fibroblast-like blastema #8","Fibroblast-like blastema #9",
  "Myogenic blastema","Myogenic blastemaa",
  "Endothelial","Pericyte","Schwann"
)

CT_LABELS <- c(
  "Macrophage"                  = "Macro.",
  "Recruited_Macrophage"        = "Recr. Macro.",
  "Phagocytosing_Neutrophil"    = "Phago. Neutro.",
  "Neutrophil_#1"               = "Neutro. #1",
  "Neutrophil_#2"               = "Neutro. #2",
  "T_cell"                      = "T cell",
  "Early_B"                     = "Early B",
  "Erythrocyte"                 = "Erythro.",
  "Erythrocyte #1"              = "Erythro. #1",
  "Erythrocyte #2"              = "Erythro. #2",
  "Intermediate_epidermis"      = "Int. Epid.",
  "Proliferating_Epidermis"     = "Prolif. Epid.",
  "BasalEpidermis"              = "Basal Epid.",
  "Epidermal_Langerhans"        = "Langerhans",
  "Basal wound epidermis"       = "Basal WE",
  "Basal WE "                   = "Basal WE",
  "Intermediate_WE_#1"          = "Int. WE #1",
  "Intermediate_WE_#2"          = "Int. WE #2",
  "Intermediate_WE_#3"          = "Int. WE #3",
  "Intermediate WE #4"          = "Int. WE #4",
  "Intermediate WE #5"          = "Int. WE #5",
  "SSC"                         = "SSC",
  "SSCs"                        = "SSC",
  "Fibroblast"                  = "Fibro.",
  "Fibroblast-like_blastema_#1" = "Blast. #1",
  "Fibroblast-like_blastema_#2" = "Blast. #2",
  "Fibroblast-like blastema #3" = "Blast. #3",
  "Fibroblast-like blastema #4" = "Blast. #4",
  "Fibroblast-like blastema #5" = "Blast. #5",
  "Fibroblast-like blastema #7" = "Blast. #7",
  "Fibroblast-like blastema #8" = "Blast. #8",
  "Fibroblast-like blastema #9" = "Blast. #9",
  "Myogenic blastema"           = "Myog. Blast.",
  "Myogenic blastemaa"          = "Myog. Blast.",
  "Endothelial"                 = "Endo.",
  "Pericyte"                    = "Peri.",
  "Schwann"                     = "Schwann"
)

GOI <- c(
  "GLI1","PTCH1","PTCH2","SHH","IHH","DHH",
  "GDF5","MSX1","SALL1","GREM1",
  "FREM1","CYTL1","CLDN5","IL11","FGFBP1",
  "C3","LDLR","ALOX15B","MSMO1","CYP51A1","FDPS"
)

# =============================================================================
# LOAD & PREPARE
# =============================================================================

message("Loading axo_integrated.rds...")
axo <- readRDS(file.path(results_dir, "axo_integrated.rds"))
axo <- JoinLayers(axo)

axo$timepoint_pooled <- as.character(axo$timepoint)
axo$timepoint_pooled[axo$timepoint_pooled %in% WH_REPS] <- "WH"
axo$timepoint_pooled <- factor(axo$timepoint_pooled, levels = TIMEPOINT_ORDER)

goi_in_sc <- intersect(GOI, rownames(axo))
message("GOI in SC object: ", length(goi_in_sc), "/", length(GOI))

# =============================================================================
# BUILD PER-(GENE x CELL_TYPE x TIMEPOINT) SUMMARY
# =============================================================================

message("Extracting expression matrix...")
mat    <- GetAssayData(axo, layer = "data")
ct_vec <- as.character(axo$paper_cluster)
tp_vec <- as.character(axo$timepoint_pooled)

ct_present <- CT_ORDER[CT_ORDER %in% unique(ct_vec)]
ct_levels  <- c(ct_present, setdiff(unique(ct_vec), ct_present))

message("Building summary (", length(ct_levels),
        " cell types x ", length(TIMEPOINT_ORDER), " timepoints)...")

sc_summary <- bind_rows(lapply(goi_in_sc, function(g) {
  expr <- as.numeric(mat[g, ])
  bind_rows(lapply(ct_levels, function(ct) {
    bind_rows(lapply(TIMEPOINT_ORDER, function(tp) {
      idx <- which(ct_vec == ct & tp_vec == tp)
      if (length(idx) == 0) return(NULL)
      e <- expr[idx]
      data.frame(
        gene      = g,
        cell_type = ct,
        timepoint = tp,
        mean_expr = mean(e),
        pct_expr  = mean(e > 0) * 100,
        stringsAsFactors = FALSE
      )
    }))
  }))
}))

sc_summary$cell_type <- factor(sc_summary$cell_type, levels = ct_levels)
sc_summary$timepoint <- factor(sc_summary$timepoint, levels = TIMEPOINT_ORDER)
sc_summary$ct_label  <- CT_LABELS[as.character(sc_summary$cell_type)]
sc_summary$ct_label[is.na(sc_summary$ct_label)] <-
  as.character(sc_summary$cell_type)[is.na(sc_summary$ct_label)]

message("Summary: ", nrow(sc_summary), " rows")

# =============================================================================
# COLOR PALETTE
# 36 biologically ordered cell types — grouped by broad lineage so adjacent
# colors in the legend reflect biological relatedness.
# =============================================================================

LINEAGE_COLORS <- c(
  # Immune (blue-green family)
  "Macrophage"                  = "#1A5276",
  "Recruited_Macrophage"        = "#1F618D",
  "Phagocytosing_Neutrophil"    = "#2E86C1",
  "Neutrophil_#1"               = "#5DADE2",
  "Neutrophil_#2"               = "#85C1E9",
  "T_cell"                      = "#AED6F1",
  "Early_B"                     = "#D6EAF8",
  "Erythrocyte"                 = "#7B241C",
  "Erythrocyte #1"              = "#A93226",
  "Erythrocyte #2"              = "#CD6155",
  # Epithelial (green family)
  "Intermediate_epidermis"      = "#1E8449",
  "Proliferating_Epidermis"     = "#27AE60",
  "BasalEpidermis"              = "#52BE80",
  "Epidermal_Langerhans"        = "#82E0AA",
  "Basal wound epidermis"       = "#145A32",
  "Basal WE "                   = "#1D6A39",
  "Intermediate_WE_#1"          = "#117A65",
  "Intermediate_WE_#2"          = "#148F77",
  "Intermediate_WE_#3"          = "#17A589",
  "Intermediate WE #4"          = "#1ABC9C",
  "Intermediate WE #5"          = "#76D7C4",
  # Fibroblast / blastema (orange-red family)
  "SSC"                         = "#784212",
  "SSCs"                        = "#784212",
  "Fibroblast"                  = "#935116",
  "Fibroblast-like_blastema_#1" = "#BA4A00",
  "Fibroblast-like_blastema_#2" = "#CB4335",
  "Fibroblast-like blastema #3" = "#E74C3C",
  "Fibroblast-like blastema #4" = "#EC7063",
  "Fibroblast-like blastema #5" = "#F1948A",
  "Fibroblast-like blastema #7" = "#D35400",
  "Fibroblast-like blastema #8" = "#E67E22",
  "Fibroblast-like blastema #9" = "#F0B27A",
  "Myogenic blastema"           = "#F4D03F",
  "Myogenic blastemaa"          = "#F4D03F",
  # Vascular / neural (purple family)
  "Endothelial"                 = "#6C3483",
  "Pericyte"                    = "#8E44AD",
  "Schwann"                     = "#BB8FCE"
)

# =============================================================================
# SHARED THEME
# =============================================================================

base_theme <- theme_classic(base_size = 10) +
  theme(
    axis.text.x        = element_text(size = 9, color = "grey20"),
    axis.text.y        = element_text(size = 8),
    axis.title         = element_text(size = 9),
    plot.title         = element_text(size = 10, face = "bold.italic",
                                      hjust = 0, margin = margin(b = 3)),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    plot.margin        = margin(t = 5, r = 5, b = 10, l = 5)
  )

leg_theme <- theme(
  legend.position  = "right",
  legend.title     = element_text(size = 7.5, face = "bold"),
  legend.text      = element_text(size = 7),
  legend.key.size  = unit(0.4, "cm")
)

# =============================================================================
# PER-GENE TEMPORAL TRAJECTORY FIGURE
# =============================================================================

message("\nGenerating temporal expression figures for ", length(goi_in_sc), " genes...")

for (g in goi_in_sc) {
  df <- sc_summary |> filter(gene == g)
  if (nrow(df) == 0) next

  # Show only cell types with detectable expression in at least one stage
  expressing_cts <- df |>
    group_by(cell_type) |>
    summarise(max_expr = max(mean_expr), .groups = "drop") |>
    filter(max_expr > 0.05) |>
    pull(cell_type)

  df_plot <- df |> filter(cell_type %in% expressing_cts)

  if (nrow(df_plot) == 0) {
    message("  Skipping ", g, " — no cell types exceed expression threshold")
    next
  }

  ct_colors_use <- LINEAGE_COLORS[as.character(unique(df_plot$cell_type))]
  ct_colors_use[is.na(ct_colors_use)] <- "grey50"

  # Left panel: mean log-normalized expression across stages
  px <- ggplot(df_plot, aes(x = timepoint, y = mean_expr,
                             color = cell_type, group = cell_type)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    geom_point(size = 2.2, alpha = 0.95) +
    scale_color_manual(
      values = ct_colors_use,
      name   = "Cell type",
      labels = function(x) ifelse(!is.na(CT_LABELS[x]), CT_LABELS[x], x)
    ) +
    scale_x_discrete(limits = TIMEPOINT_ORDER) +
    labs(title = g,
         x     = "Regeneration stage",
         y     = "Mean log-norm. expression") +
    base_theme + leg_theme

  # Right panel: percent of cells expressing
  pp <- ggplot(df_plot, aes(x = timepoint, y = pct_expr,
                             color = cell_type, group = cell_type)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    geom_point(size = 2.2, alpha = 0.95) +
    scale_color_manual(
      values = ct_colors_use,
      name   = "Cell type",
      labels = function(x) ifelse(!is.na(CT_LABELS[x]), CT_LABELS[x], x)
    ) +
    scale_x_discrete(limits = TIMEPOINT_ORDER) +
    labs(title = paste0(g, " - % expressing"),
         x     = "Regeneration stage",
         y     = "% cells expressing") +
    base_theme + leg_theme

  n_ct <- length(expressing_cts)
  page <- px + pp + plot_layout(ncol = 2, widths = c(1, 1))

  pdf(file.path(out_dir, paste0(g, "_expr_curve.pdf")),
      width = 14, height = max(4, 3 + n_ct * 0.1))
  print(page)
  dev.off()
  message("  Saved: ", g, " (", n_ct, " expressing cell types)")
}

message("\nAll expression curve figures saved to ", out_dir)
