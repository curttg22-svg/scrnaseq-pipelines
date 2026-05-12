# =============================================================================
# plot_hedgehog_xen_axo_comparison.R
#
# Side-by-side violin plots comparing Ptch1, Ptch2, Gli1 expression across
# matched cell populations in:
#   - Xenopus froglet 7-14 dpa (Lin et al. 2021)  — non-regenerative
#   - Axolotl mid-blastema / MB (Leigh et al. 2018) — regenerative
#
# PTCH2 is absent from the axolotl dataset (not detected).
# Run from comparison/ directory: source("R/plot_hedgehog_xen_axo_comparison.R")
#
# Output: results/hedgehog_xen_axo_comparison.pdf
# =============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)

XEN_RDS   <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
AXO_RDS   <- "~/Desktop/Axolotl_scrnaseq_project_V3/results/axo_integrated.rds"
RESULTS_DIR <- "results"

# =============================================================================
# HELPERS
# =============================================================================

sum_homeologs <- function(expr_mat, pattern) {
  hits <- grep(pattern, rownames(expr_mat), value = TRUE, perl = TRUE, ignore.case = TRUE)
  if (length(hits) == 0) return(rep(0, ncol(expr_mat)))
  if (length(hits) == 1) return(as.numeric(expr_mat[hits, ]))
  Matrix::colSums(expr_mat[hits, , drop = FALSE])
}

CELL_COLORS <- c(
  "Keratinocyte"   = "#6A1B9A",
  "Basal epidermis"= "#AB47BC",
  "Macrophage"     = "#E65100",
  "CT fibroblast"  = "#1565C0",
  "CT progenitor"  = "#42A5F5"
)

violin_theme <- theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 9),
    strip.text       = element_text(face = "bold.italic", size = 11),
    strip.background = element_rect(fill = "grey93", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, color = "grey40")
  )

make_violin <- function(df, genes, title, subtitle) {
  df_long <- data.frame(
    cell_label = rep(df$cell_label, length(genes)),
    gene       = rep(genes, each = nrow(df)),
    expression = unlist(lapply(genes, function(g) df[[g]]))
  )
  df_long <- df_long[!is.na(df_long$cell_label), ]
  df_long$gene <- factor(df_long$gene, levels = genes)

  ggplot(df_long, aes(x = cell_label, y = expression, fill = cell_label)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.35, color = "grey30") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 color = "grey20", linewidth = 0.4) +
    scale_fill_manual(values = CELL_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    facet_wrap(~ gene, nrow = 1, scales = "free_y") +
    labs(title = title, subtitle = subtitle, x = NULL,
         y = "Expression (summed log-norm)") +
    violin_theme
}

# =============================================================================
# 1. XENOPUS — LOAD ONCE, JOIN LAYERS, EXTRACT BOTH CONDITIONS
# =============================================================================

XEN_CT <- c("Keratinocyte", "Epithelial_mixed",
            "Macrophage", "Macrophage_2",
            "CT_fibroblast_A", "CT_fibroblast_B", "CT_fibroblast_C",
            "CT_fibroblast_D", "CT_fibroblast_E", "CT_progenitor")

label_xen <- function(cell_type_vec) {
  grp <- as.character(cell_type_vec)
  grp[grepl("^CT_fibroblast", grp)] <- "CT_fibroblast"
  grp[grp == "Macrophage_2"]        <- "Macrophage"
  factor(grp,
         levels = c("Keratinocyte","Epithelial_mixed","Macrophage",
                    "CT_fibroblast","CT_progenitor"),
         labels = c("Keratinocyte","Basal epidermis","Macrophage",
                    "CT fibroblast","CT progenitor"))
}

xen_extract <- function(rds_path, cond, tp, ct_types) {
  xen <- readRDS(path.expand(rds_path))
  # Subset by metadata first — works on split layers, much smaller object
  sub <- subset(xen, condition == cond & timepoint == tp & cell_type %in% ct_types)
  rm(xen); gc()
  # Join layers only on the small subset
  sub <- JoinLayers(sub)
  expr <- GetAssayData(sub, layer = "data")
  df <- data.frame(
    cell_label = label_xen(sub$cell_type),
    PTCH1 = sum_homeologs(expr, "(?i)^ptch1[.]"),
    PTCH2 = sum_homeologs(expr, "(?i)^ptch2[.]"),
    GLI1  = sum_homeologs(expr, "(?i)^gli1[.]")
  )
  df[!is.na(df$cell_label), ]
}

message("Loading Xenopus froglet 7-14 dpa...")
xen_df  <- xen_extract(XEN_RDS, "BL",   "7-14dpa", XEN_CT)

message("Loading Xenopus tadpole NF52...")
lbst_df <- xen_extract(XEN_RDS, "LBst", "NF52", XEN_CT)

# Pre-load each tadpole timepoint for individual PDFs
message("Loading Xenopus tadpole NF50...")
lbst_NF50 <- xen_extract(XEN_RDS, "LBst", "NF50", XEN_CT)
message("Loading Xenopus tadpole NF51...")
lbst_NF51 <- xen_extract(XEN_RDS, "LBst", "NF51", XEN_CT)
lbst_NF52 <- lbst_df  # already loaded

xen_n <- table(xen_df$cell_label)
xen_subtitle <- paste0("n = ", paste(paste0(names(xen_n), " (", xen_n, ")"),
                                     collapse = ", "))

p_xen <- make_violin(xen_df, c("PTCH1","PTCH2","GLI1"),
                     title    = "Xenopus froglet 7-14 dpa (Lin et al. 2021)",
                     subtitle = xen_subtitle)

# =============================================================================
# 2. AXOLOTL — EARLY AND MID BLASTEMA
# =============================================================================

message("Loading Axolotl...")
axo <- readRDS(path.expand(AXO_RDS))
axo$paper_cluster <- trimws(axo$paper_cluster)

AXO_CLUSTERS <- c("Basal WE", "Intermediate WE", "Macrophage",
                  "Fibroblast-like blastema")
AXO_TP_ORDER <- c("EB", "MB")
AXO_TP_LABELS <- c("EB" = "Early blastema", "MB" = "Mid-blastema")

make_axo_panel <- function(axo, tp) {
  sub <- subset(axo, timepoint == tp & paper_cluster %in% AXO_CLUSTERS)
  expr <- GetAssayData(sub, layer = "data")
  df <- data.frame(
    cell_label = factor(
      trimws(sub$paper_cluster),
      levels = c("Intermediate WE", "Basal WE", "Macrophage",
                 "Fibroblast-like blastema"),
      labels = c("Keratinocyte", "Basal epidermis", "Macrophage", "CT progenitor")
    ),
    PTCH1 = as.numeric(expr["PTCH1", ]),
    GLI1  = as.numeric(expr["GLI1",  ])
  )
  df <- df[!is.na(df$cell_label), ]
  n   <- table(df$cell_label)
  sub_title <- paste0("n = ", paste(paste0(names(n), " (", n, ")"), collapse = ", "),
                      "  |  PTCH2 not detected")
  make_violin(df, c("PTCH1", "GLI1"),
              title    = paste0("Axolotl ", AXO_TP_LABELS[tp], " (Leigh et al. 2018)"),
              subtitle = sub_title)
}

p_mb <- make_axo_panel(axo, "MB")
rm(axo); gc()

# (LBst data extracted above alongside BL)

lbst_n <- table(lbst_df$cell_label)
lbst_subtitle <- paste0("n = ", paste(paste0(names(lbst_n), " (", lbst_n, ")"),
                                       collapse = ", "),
                         "  |  * macrophage n low")

p_lbst <- make_violin(lbst_df, c("PTCH1", "PTCH2", "GLI1"),
                      title    = "Regenerative limb bud",
                      subtitle = paste0("Xenopus tadpole NF52 (Lin et al. 2021)  |  ",
                                        lbst_subtitle))

# =============================================================================
# 3. COMBINE AND SAVE
# =============================================================================

p_xen_v2 <- p_xen +
  labs(title    = "Non-regenerative wound healing",
       subtitle = "Xenopus froglet 7-14 dpa (Lin et al. 2021)")

p_mb_v2 <- p_mb +
  labs(title    = "Regenerative blastema",
       subtitle = "Axolotl mid-blastema (Leigh et al. 2018)  |  PTCH2 not detected")

combined <- p_xen_v2 / p_lbst / p_mb_v2 +
  plot_annotation(
    title    = "Hedgehog pathway gene expression across regenerative contexts",
    subtitle = "Homeologs summed per cell | log-normalized expression",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  )

# Three-panel version (froglet + tadpole + axolotl MB)
out_file <- file.path(RESULTS_DIR, "hedgehog_xen_axo_MB_only.pdf")
pdf(out_file, width = 11, height = 13)
print(combined)
dev.off()
message("Saved ", out_file)

# Two-panel version (froglet + axolotl MB only — no tadpole)
combined_2panel <- p_xen_v2 / p_mb_v2 +
  plot_annotation(
    title    = "Hedgehog pathway gene expression across regenerative contexts",
    subtitle = "Homeologs summed per cell | log-normalized expression",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  )

out_file2 <- file.path(RESULTS_DIR, "hedgehog_nonregen_vs_blastema.pdf")
pdf(out_file2, width = 11, height = 9)
print(combined_2panel)
dev.off()
message("Saved ", out_file2)

# =============================================================================
# Individual PDFs: froglet 7-14dpa vs each tadpole timepoint
# =============================================================================

tp_data <- list(NF50 = lbst_NF50, NF51 = lbst_NF51, NF52 = lbst_NF52)

for (tp in names(tp_data)) {
  tp_n  <- table(tp_data[[tp]]$cell_label)
  tp_sub <- paste0("n = ", paste(paste0(names(tp_n), " (", tp_n, ")"),
                                  collapse = ", "),
                   "  |  * macrophage n low")

  p_tp <- make_violin(tp_data[[tp]], c("PTCH1", "PTCH2", "GLI1"),
                      title    = "Regenerative limb bud",
                      subtitle = paste0("Xenopus tadpole ", tp,
                                        " (Lin et al. 2021)  |  ", tp_sub))

  combo <- p_xen_v2 / p_tp +
    plot_annotation(
      title    = paste0("Hedgehog pathway: non-regenerative wound healing vs tadpole ", tp),
      subtitle = "Homeologs summed per cell | log-normalized expression",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40")
      )
    )

  fn <- file.path(RESULTS_DIR, paste0("hedgehog_froglet_7-14dpa_vs_", tp, ".pdf"))
  pdf(fn, width = 11, height = 9)
  print(combo)
  dev.off()
  message("Saved ", fn)
}
