# =============================================================================
# 04 — LOESS expression curves: mean expression per cell type across timepoints
# Input:  results/axo_integrated.rds
# Output: results/loess_curves/<GENE>_loess.pdf   (one per GOI)
#
# WH_N4 / WH_N5 / WH_N6 are pooled into "WH" before computing means,
# so the final x-axis shows four stages: Intact → WH → EB → MB.
# Right panel: percent of cells expressing the gene (dot size proxy).
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
out_dir     <- file.path(results_dir, "loess_curves")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

TIMEPOINT_ORDER <- c("Intact", "WH", "EB", "MB")
WH_REPS         <- c("WH_N4", "WH_N5", "WH_N6")

TIMEPOINT_COLORS <- c(
  "Intact" = "#2E7D32",
  "WH"     = "#00838F",
  "EB"     = "#1565C0",
  "MB"     = "#6A1B9A"
)

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

# Pool WH replicates into a single "WH" stage
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
sc_summary$ct_num <- as.integer(sc_summary$cell_type)

message("Summary: ", nrow(sc_summary), " rows")

# =============================================================================
# SHARED THEME
# =============================================================================

base_theme <- theme_classic(base_size = 10) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1,
                                      size = 6.5, color = "grey20"),
    axis.text.y        = element_text(size = 8),
    axis.title.y       = element_text(size = 8.5),
    axis.title.x       = element_blank(),
    plot.title         = element_text(size = 10, face = "bold.italic",
                                      hjust = 0, margin = margin(b = 3)),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
    legend.position    = "none",
    plot.margin        = margin(t = 5, r = 5, b = 20, l = 5)
  )

leg_theme <- theme(
  legend.position  = "right",
  legend.title     = element_text(size = 8, face = "bold"),
  legend.text      = element_text(size = 8),
  legend.key.size  = unit(0.45, "cm")
)

# =============================================================================
# PER-GENE LOESS FIGURE
# =============================================================================

message("\nGenerating LOESS figures for ", length(goi_in_sc), " genes...")

for (g in goi_in_sc) {
  df <- sc_summary |> filter(gene == g)
  if (nrow(df) == 0) next

  ct_label_vec <- df |>
    distinct(ct_num, ct_label) |>
    arrange(ct_num)

  # Left panel: mean expression per cell type per timepoint
  # Connected dot plot — cell type ordering is categorical, not continuous,
  # so LOESS smoothing over integer positions fabricates trends between
  # unrelated cell types and is not appropriate here.
  px <- ggplot(df, aes(x = ct_num, y = mean_expr,
                        color = timepoint, group = timepoint)) +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(size = 1.8, alpha = 0.85) +
    scale_color_manual(values = TIMEPOINT_COLORS,
                       name  = "Timepoint",
                       labels = TIMEPOINT_ORDER) +
    scale_x_continuous(
      breaks = ct_label_vec$ct_num,
      labels = ct_label_vec$ct_label,
      expand = expansion(mult = 0.02)
    ) +
    labs(title = g,
         y     = "Mean log-norm. expr.") +
    base_theme + leg_theme

  # Right panel: percent expressing (connected dot plot, same rationale)
  pp <- ggplot(df, aes(x = ct_num, y = pct_expr,
                        color = timepoint, group = timepoint)) +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(size = 1.8, alpha = 0.85) +
    scale_color_manual(values = TIMEPOINT_COLORS,
                       name  = "Timepoint",
                       labels = TIMEPOINT_ORDER) +
    scale_x_continuous(
      breaks = ct_label_vec$ct_num,
      labels = ct_label_vec$ct_label,
      expand = expansion(mult = 0.02)
    ) +
    labs(title = paste0(g, " — % expressing"),
         y     = "% cells expressing") +
    base_theme + leg_theme

  page <- px + pp + plot_layout(ncol = 2, widths = c(1, 1))

  pdf(file.path(out_dir, paste0(g, "_loess.pdf")), width = 16, height = 4.2)
  print(page)
  dev.off()
  message("  Saved: ", g)
}

message("\nAll LOESS figures saved to ", out_dir)
