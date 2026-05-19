# =============================================================================
# plot_pca_exploration.R
#
# Produces individual PCA scatter plots (one file per PC pair) for the
# pseudobulk reference samples only — no BMS/EtOH data.
#
# Two approaches:
#   (A) Spearman MDS (k=5) — the same basis as the existing spectrum plots
#   (B) Direct PCA on the pseudobulk expression matrix
#
# Three normalizations: raw / within-species centered / within-species z-scored
#
# Each PC pair saved as its own PNG in figures/pca_exploration/.
# Also saves a scree plot per approach × normalization.
#
# Run from comparison/ directory: source("R/plot_pca_exploration.R")
# =============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

CACHE   <- "results/pseudobulk_cache.rds"
FIG_DIR <- "figures/pca_exploration"
dir.create(FIG_DIR, showWarnings = FALSE)

# =============================================================================
# 1. Load cache
# =============================================================================

cache    <- readRDS(CACHE)
common   <- cache$common_genes
pseudos  <- cache$pseudos

axo_full <- pseudos[["Axolotl (collaborator)"]]
mus_full <- pseudos[["Mouse digit tip"]]
xen_all  <- pseudos[["Xenopus (Lin 2021)"]]

axo_mat  <- axo_full[common, ]
mus_mat  <- mus_full[common, ]
xen_keep <- !colnames(xen_all) %in% c("LinBL_7-14dpa", "LinBL_14-52dpa")
xen_mat  <- xen_all[common, xen_keep, drop = FALSE]

ref_mat <- cbind(axo_mat, xen_mat, mus_mat)

ref_sample_species <- c(
  setNames(rep("Axolotl", ncol(axo_mat)), colnames(axo_mat)),
  setNames(rep("Xenopus", ncol(xen_mat)), colnames(xen_mat)),
  setNames(rep("Mouse",   ncol(mus_mat)), colnames(mus_mat))
)

# =============================================================================
# 2. Normalisation
# =============================================================================

normalize_ref <- function(mat, sp_vec, method) {
  if (method == "none") return(mat)
  result <- mat
  for (sp in unique(sp_vec)) {
    cols <- names(sp_vec)[sp_vec == sp]
    m    <- rowMeans(mat[, cols, drop = FALSE])
    if (method == "center") {
      result[, cols] <- mat[, cols] - m
    } else {
      s <- apply(mat[, cols, drop = FALSE], 1, sd); s[s < 1e-9] <- 1
      result[, cols] <- (mat[, cols] - m) / s
    }
  }
  result
}

# =============================================================================
# 3. Sample metadata
# =============================================================================

make_group <- function(nms) {
  dplyr::case_when(
    grepl("NonRegen", nms) & grepl("Axo", nms) ~ "Axolotl (non-regen)",
    grepl("Regen",    nms) & grepl("Axo", nms) ~ "Axolotl (regenerating)",
    grepl("^LinBL_",  nms)                      ~ "Xenopus froglet (non-regen)",
    grepl("^LinLBst_",nms)                      ~ "Xenopus tadpole (regen)",
    grepl("^Mouse_",  nms)                      ~ "Mouse digit tip",
    TRUE                                         ~ "Unknown"
  )
}

make_label <- function(nms) {
  nms |>
    gsub("AxoCollab_", "",       x = _) |>
    gsub("LinBL_",     "BL ",    x = _) |>
    gsub("LinLBst_",   "LBst ",  x = _) |>
    gsub("Mouse_",     "Mouse ", x = _) |>
    gsub("_",          " ",      x = _)
}

group_colors <- c(
  "Axolotl (regenerating)"      = "#1565C0",
  "Axolotl (non-regen)"         = "#90CAF9",
  "Xenopus tadpole (regen)"     = "#00897B",
  "Xenopus froglet (non-regen)" = "#F4511E",
  "Mouse digit tip"             = "#43A047"
)
group_shapes <- c(
  "Axolotl (regenerating)"      = 21,
  "Axolotl (non-regen)"         = 24,
  "Xenopus tadpole (regen)"     = 21,
  "Xenopus froglet (non-regen)" = 21,
  "Mouse digit tip"             = 21
)

# =============================================================================
# 4. Plot helpers
# =============================================================================

make_scatter <- function(df, xcol, ycol, pct_x, pct_y,
                          title_str, subtitle_str) {
  ggplot(df, aes(.data[[xcol]], .data[[ycol]],
                  fill = group, shape = group, label = label)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
    stat_ellipse(aes(color = group, group = group),
                 type = "norm", level = 0.8,
                 linewidth = 0.35, linetype = "dashed", alpha = 0.5,
                 show.legend = FALSE) +
    geom_point(size = 5, color = "white", stroke = 0.5, alpha = 0.95) +
    geom_text_repel(size = 3, color = "grey20", bg.color = "white", bg.r = 0.08,
                    box.padding = 0.5, point.padding = 0.3,
                    max.overlaps = 25,
                    segment.color = "grey60", segment.size = 0.3) +
    scale_fill_manual(values  = group_colors, name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    scale_shape_manual(values = group_shapes, name = NULL) +
    labs(
      title    = title_str,
      subtitle = subtitle_str,
      x = sprintf("%s  (%.1f%%)", xcol, pct_x),
      y = sprintf("%s  (%.1f%%)", ycol, pct_y)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 14, color = "grey15"),
      plot.subtitle    = element_text(size = 9,  color = "grey45"),
      axis.title       = element_text(size = 11, color = "grey25"),
      axis.text        = element_text(size = 9,  color = "grey40"),
      legend.text      = element_text(size = 10, color = "grey20"),
      legend.key       = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.position  = "right",
      plot.margin      = margin(16, 16, 16, 16)
    ) +
    guides(fill  = guide_legend(override.aes = list(size = 4.5)),
           shape = guide_legend(override.aes = list(size = 4.5)))
}

make_scree <- function(pct_var, n_show = 8, title_str, subtitle_str) {
  df <- data.frame(
    PC  = factor(paste0("PC", seq_len(n_show)), levels = paste0("PC", seq_len(n_show))),
    pct = pct_var[seq_len(n_show)]
  )
  ggplot(df, aes(PC, pct)) +
    geom_col(fill = "#5C6BC0", color = "white", linewidth = 0.4, width = 0.65) +
    geom_text(aes(label = sprintf("%.1f%%", pct)),
              vjust = -0.4, size = 3.2, color = "grey30") +
    labs(title = title_str, subtitle = subtitle_str,
         x = NULL, y = "% variance explained") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
      plot.title    = element_text(face = "bold", size = 13, color = "grey15"),
      plot.subtitle = element_text(size = 9, color = "grey45"),
      axis.text     = element_text(size = 10, color = "grey35"),
      plot.margin   = margin(16, 16, 16, 16)
    )
}

# Pairs to produce individual plots for
pc_pairs <- list(c(1,2), c(1,3), c(2,3), c(1,4), c(2,4), c(3,4))

n_dims <- 6   # number of PCs / MDS dimensions to compute

# =============================================================================
# 5. Main loop
# =============================================================================

for (method in c("none", "center", "zscore")) {
  suffix     <- switch(method, none = "raw", center = "centered", zscore = "zscored")
  norm_label <- switch(method,
                       none   = "Raw",
                       center = "Centered",
                       zscore = "Z-scored")

  ref_n <- normalize_ref(ref_mat, ref_sample_species, method)

  # ---- A. Spearman MDS -------------------------------------------------------
  message("Spearman MDS | ", norm_label)
  rho_ref <- cor(ref_n, method = "spearman")
  mds_fit <- cmdscale(as.dist(1 - rho_ref), k = n_dims, eig = TRUE)
  eig_pos <- pmax(mds_fit$eig, 0)
  pct_mds <- 100 * eig_pos[seq_len(n_dims)] / sum(eig_pos)

  mds_df <- as.data.frame(mds_fit$points)
  colnames(mds_df) <- paste0("PC", seq_len(n_dims))
  mds_df$sample <- rownames(mds_fit$points)
  mds_df$label  <- make_label(mds_df$sample)
  mds_df$group  <- make_group(mds_df$sample)

  # Scree
  p_scree <- make_scree(pct_mds,
                         title_str    = sprintf("Scree — Spearman MDS | %s", norm_label),
                         subtitle_str = sprintf("%d common genes × %d samples",
                                                length(common), ncol(ref_n)))
  ggsave(file.path(FIG_DIR, paste0("mds_scree_", suffix, ".png")),
         p_scree, width = 8, height = 5, dpi = 200, bg = "white")

  # Individual PC pair plots
  for (pair in pc_pairs) {
    i <- pair[1]; j <- pair[2]
    xcol <- paste0("PC", i); ycol <- paste0("PC", j)
    p <- make_scatter(
      mds_df, xcol, ycol, pct_mds[i], pct_mds[j],
      title_str    = sprintf("Spearman MDS — %s vs %s | %s", xcol, ycol, norm_label),
      subtitle_str = sprintf("Spearman distance MDS | %d common genes | %s normalisation",
                             length(common), norm_label)
    )
    fname <- file.path(FIG_DIR,
                       sprintf("mds_%s_PC%dvPC%d.png", suffix, i, j))
    ggsave(fname, p, width = 9, height = 7, dpi = 200, bg = "white")
    message("  Saved: ", basename(fname))
  }

  # ---- B. Direct PCA ---------------------------------------------------------
  message("Direct PCA   | ", norm_label)
  pca_fit <- prcomp(t(ref_n), center = TRUE, scale. = FALSE)
  pct_pca <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)

  pca_df <- as.data.frame(pca_fit$x[, seq_len(n_dims)])
  colnames(pca_df) <- paste0("PC", seq_len(n_dims))
  pca_df$sample <- rownames(pca_fit$x)
  pca_df$label  <- make_label(pca_df$sample)
  pca_df$group  <- make_group(pca_df$sample)

  # Scree
  p_scree_pca <- make_scree(pct_pca,
                              title_str    = sprintf("Scree — Direct PCA | %s", norm_label),
                              subtitle_str = sprintf("%d common genes × %d samples",
                                                     length(common), ncol(ref_n)))
  ggsave(file.path(FIG_DIR, paste0("pca_scree_", suffix, ".png")),
         p_scree_pca, width = 8, height = 5, dpi = 200, bg = "white")

  # Individual PC pair plots
  for (pair in pc_pairs) {
    i <- pair[1]; j <- pair[2]
    xcol <- paste0("PC", i); ycol <- paste0("PC", j)
    p <- make_scatter(
      pca_df, xcol, ycol, pct_pca[i], pct_pca[j],
      title_str    = sprintf("Direct PCA — %s vs %s | %s", xcol, ycol, norm_label),
      subtitle_str = sprintf("PCA on pseudobulk expression | %d common genes | %s normalisation",
                             length(common), norm_label)
    )
    fname <- file.path(FIG_DIR,
                       sprintf("pca_%s_PC%dvPC%d.png", suffix, i, j))
    ggsave(fname, p, width = 9, height = 7, dpi = 200, bg = "white")
    message("  Saved: ", basename(fname))
  }
}

message("\nDone. ", length(list.files(FIG_DIR, "*.png")),
        " figures saved to ", FIG_DIR, "/")
