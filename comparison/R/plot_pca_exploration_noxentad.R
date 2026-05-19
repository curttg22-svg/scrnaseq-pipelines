# =============================================================================
# plot_pca_exploration_noxentad.R
#
# Same as plot_pca_exploration.R but Xenopus tadpole (LinLBst_) samples are
# excluded. Reference set = Axolotl + Xenopus froglet (LinBL_) + Mouse only.
#
# Also produces PC1vsPC2, PC2vsPC3, PC2vsPC4 with BMS / EtOH projected.
#
# Run from comparison/ directory: source("R/plot_pca_exploration_noxentad.R")
# =============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

CACHE   <- "results/pseudobulk_cache.rds"
FIG_DIR <- "figures/pca_noxentad"
dir.create(FIG_DIR, showWarnings = FALSE)

# =============================================================================
# 1. Load cache — drop Xenopus tadpole (LinLBst_) samples
# =============================================================================

cache    <- readRDS(CACHE)
common   <- cache$common_genes
pseudos  <- cache$pseudos
bulk_mat <- cache$bulk_mat

axo_full <- pseudos[["Axolotl (collaborator)"]]
mus_full <- pseudos[["Mouse digit tip"]]
xen_all  <- pseudos[["Xenopus (Lin 2021)"]]

axo_mat <- axo_full[common, ]
mus_mat <- mus_full[common, ]

# Keep only froglet (LinBL_) Xenopus; drop tadpole (LinLBst_) and the
# two aggregate timepoint bins excluded previously
xen_keep <- grepl("^LinBL_", colnames(xen_all)) &
            !colnames(xen_all) %in% c("LinBL_7-14dpa", "LinBL_14-52dpa")
xen_mat  <- xen_all[common, xen_keep, drop = FALSE]
xen_full <- xen_all[, xen_keep, drop = FALSE]

ref_mat <- cbind(axo_mat, xen_mat, mus_mat)

ref_sample_species <- c(
  setNames(rep("Axolotl", ncol(axo_mat)), colnames(axo_mat)),
  setNames(rep("Xenopus", ncol(xen_mat)), colnames(xen_mat)),
  setNames(rep("Mouse",   ncol(mus_mat)), colnames(mus_mat))
)
ref_species_mats <- list(Axolotl = axo_full, Xenopus = xen_full, Mouse = mus_full)

cat(sprintf("Reference samples: Axolotl=%d  Xenopus froglet=%d  Mouse=%d  (total=%d)\n",
            ncol(axo_mat), ncol(xen_mat), ncol(mus_mat), ncol(ref_mat)))

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
    grepl("^Mouse_",  nms)                      ~ "Mouse digit tip",
    TRUE                                         ~ "Unknown"
  )
}

make_label <- function(nms) {
  nms |>
    gsub("AxoCollab_", "",       x = _) |>
    gsub("LinBL_",     "BL ",    x = _) |>
    gsub("Mouse_",     "Mouse ", x = _) |>
    gsub("_",          " ",      x = _)
}

group_colors <- c(
  "Axolotl (regenerating)"      = "#1565C0",
  "Axolotl (non-regen)"         = "#90CAF9",
  "Xenopus froglet (non-regen)" = "#F4511E",
  "Mouse digit tip"             = "#43A047",
  "BMS (axolotl limb)"          = "#E65100",
  "EtOH (axolotl limb)"         = "#6A1B9A"
)
group_shapes <- c(
  "Axolotl (regenerating)"      = 21,
  "Axolotl (non-regen)"         = 24,
  "Xenopus froglet (non-regen)" = 21,
  "Mouse digit tip"             = 21,
  "BMS (axolotl limb)"          = 23,
  "EtOH (axolotl limb)"         = 22
)

# =============================================================================
# 4. Plot helpers
# =============================================================================

make_scatter <- function(df, xcol, ycol, pct_x, pct_y, title_str, subtitle_str) {
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
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 14, color = "grey15"),
      plot.subtitle     = element_text(size = 9,  color = "grey45"),
      axis.title        = element_text(size = 11, color = "grey25"),
      axis.text         = element_text(size = 9,  color = "grey40"),
      legend.text       = element_text(size = 10, color = "grey20"),
      legend.key        = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.position   = "right",
      plot.margin       = margin(16, 16, 16, 16)
    ) +
    guides(fill  = guide_legend(override.aes = list(size = 4.5)),
           shape = guide_legend(override.aes = list(size = 4.5)))
}

make_scatter_bms <- function(all_df, seg_df, ring_df,
                              xcol, ycol, pct_x, pct_y,
                              title_str, subtitle_str) {
  ggplot(all_df, aes(.data[[xcol]], .data[[ycol]],
                      fill = group, shape = group, label = label)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
    geom_segment(data = seg_df,
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 inherit.aes = FALSE,
                 color = "grey50", linewidth = 0.45, linetype = "dashed") +
    geom_point(data = ring_df, aes(.data[[xcol]], .data[[ycol]]),
               inherit.aes = FALSE,
               size = 9, shape = 21, fill = NA, color = "#FF6F00", stroke = 1.8) +
    geom_point(size = 5, color = "white", stroke = 0.5, alpha = 0.95) +
    geom_text_repel(size = 3, color = "grey20", bg.color = "white", bg.r = 0.08,
                    box.padding = 0.5, point.padding = 0.3, max.overlaps = 30,
                    segment.color = "grey60", segment.size = 0.3) +
    scale_fill_manual(values  = group_colors, name = NULL) +
    scale_shape_manual(values = group_shapes, name = NULL) +
    labs(
      title    = title_str,
      subtitle = subtitle_str,
      x = sprintf("%s  (%.1f%%)", xcol, pct_x),
      y = sprintf("%s  (%.1f%%)", ycol, pct_y)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 14, color = "grey15"),
      plot.subtitle     = element_text(size = 9,  color = "grey45"),
      axis.title        = element_text(size = 11, color = "grey25"),
      axis.text         = element_text(size = 9,  color = "grey40"),
      legend.text       = element_text(size = 10, color = "grey20"),
      legend.key        = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.position   = "right",
      plot.margin       = margin(16, 16, 16, 16)
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
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
      plot.title    = element_text(face = "bold", size = 13, color = "grey15"),
      plot.subtitle = element_text(size = 9, color = "grey45"),
      axis.text     = element_text(size = 10, color = "grey35"),
      plot.margin   = margin(16, 16, 16, 16)
    )
}

pc_pairs <- list(c(1,2), c(1,3), c(2,3), c(1,4), c(2,4), c(3,4))
n_dims   <- 6

# =============================================================================
# 5. BMS / EtOH projection helper (shared Spearman centroid approach)
# =============================================================================

project_centroid_6d <- function(new_mat, sp_c, ref_sample_species, ref_coords) {
  ref_samps <- rownames(ref_coords)
  rho <- matrix(0, nrow = length(ref_samps), ncol = ncol(new_mat),
                dimnames = list(ref_samps, colnames(new_mat)))
  for (r in ref_samps) {
    sp_mat     <- sp_c[[ref_sample_species[r]]]
    pair_genes <- intersect(rownames(sp_mat), rownames(new_mat))
    rho[r, ]   <- cor(sp_mat[pair_genes, r, drop = FALSE],
                      new_mat[pair_genes, , drop = FALSE],
                      method = "spearman")[1, ]
  }
  vapply(colnames(new_mat), function(s) {
    w <- pmax(rho[, s], 0)
    if (sum(w) == 0) return(rep(0, n_dims))
    colSums(ref_coords * w) / sum(w)
  }, numeric(n_dims))
}

bms_cols  <- grep("^BMS_|^Pretreatment$|^Reamputation$", colnames(bulk_mat), value = TRUE)
etoh_cols <- grep("^EtOH_", colnames(bulk_mat), value = TRUE)

cond_pairs <- data.frame(
  bms_samp  = c("BMS_acute","BMS_24hrs","Pretreatment","Reamputation"),
  etoh_samp = c("EtOH_acute","EtOH_24hrs","EtOH_pretreat","EtOH_reamp"),
  stringsAsFactors = FALSE
)

baseline_samps <- c("AxoCollab_NonRegen_limb","LinBL_0dpa","Mouse_0dpa")

# =============================================================================
# 6. Main loop: normalisation × approach
# =============================================================================

for (method in c("none", "center", "zscore")) {
  suffix     <- switch(method, none = "raw", center = "centered", zscore = "zscored")
  norm_label <- switch(method, none = "Raw", center = "Centered", zscore = "Z-scored")

  ref_n  <- normalize_ref(ref_mat, ref_sample_species, method)

  # Centered versions of per-species matrices for projection (only used in BMS plots)
  sp_c <- if (method == "center") {
    lapply(ref_species_mats, function(m) m - rowMeans(m))
  } else if (method == "zscore") {
    lapply(ref_species_mats, function(m) {
      mu <- rowMeans(m); s <- apply(m, 1, sd); s[s < 1e-9] <- 1
      (m - mu) / s
    })
  } else {
    ref_species_mats
  }

  # ---- A. Spearman MDS -------------------------------------------------------
  message("Spearman MDS | ", norm_label)
  rho_ref <- cor(ref_n, method = "spearman")
  mds_fit <- cmdscale(as.dist(1 - rho_ref), k = n_dims, eig = TRUE)
  eig_pos <- pmax(mds_fit$eig, 0)
  pct_mds <- 100 * eig_pos[seq_len(n_dims)] / sum(eig_pos)

  ref_coords <- mds_fit$points
  colnames(ref_coords) <- paste0("PC", seq_len(n_dims))

  mds_df <- as.data.frame(ref_coords)
  mds_df$sample <- rownames(ref_coords)
  mds_df$label  <- make_label(mds_df$sample)
  mds_df$group  <- make_group(mds_df$sample)

  # Scree
  p_scree <- make_scree(pct_mds,
                         title_str    = sprintf("Scree — Spearman MDS | %s (no tadpole)", norm_label),
                         subtitle_str = sprintf("%d common genes × %d samples (Xenopus tadpole excluded)",
                                                length(common), ncol(ref_n)))
  ggsave(file.path(FIG_DIR, paste0("mds_scree_", suffix, ".png")),
         p_scree, width = 8, height = 5, dpi = 200, bg = "white")

  # Individual PC pair plots (reference only)
  for (pair in pc_pairs) {
    i <- pair[1]; j <- pair[2]
    xcol <- paste0("PC", i); ycol <- paste0("PC", j)
    p <- make_scatter(
      mds_df, xcol, ycol, pct_mds[i], pct_mds[j],
      title_str    = sprintf("Spearman MDS — %s vs %s | %s (no tadpole)", xcol, ycol, norm_label),
      subtitle_str = sprintf("Centered Spearman MDS | %d common genes | %s | Xenopus tadpole excluded",
                             length(common), norm_label)
    )
    fname <- file.path(FIG_DIR, sprintf("mds_%s_PC%dvPC%d.png", suffix, i, j))
    ggsave(fname, p, width = 9, height = 7, dpi = 200, bg = "white")
    message("  Saved: ", basename(fname))
  }

  # BMS / EtOH projection for key PC pairs (PC1v2, PC2v3, PC2v4)
  message("  Projecting BMS/EtOH | ", norm_label)
  bms_6d  <- project_centroid_6d(bulk_mat[, bms_cols,  drop = FALSE],
                                  sp_c, ref_sample_species, ref_coords)
  etoh_6d <- project_centroid_6d(bulk_mat[, etoh_cols, drop = FALSE],
                                  sp_c, ref_sample_species, ref_coords)

  for (pair in list(c(1,2), c(2,3), c(2,4))) {
    i <- pair[1]; j <- pair[2]
    xcol <- paste0("PC", i); ycol <- paste0("PC", j)

    bms_df_p <- data.frame(
      sample = colnames(bms_6d), group = "BMS (axolotl limb)",
      label  = gsub("BMS_","BMS ",gsub("_"," ",colnames(bms_6d))),
      stringsAsFactors = FALSE
    )
    bms_df_p[[xcol]] <- bms_6d[i, ]
    bms_df_p[[ycol]] <- bms_6d[j, ]

    etoh_df_p <- data.frame(
      sample = colnames(etoh_6d), group = "EtOH (axolotl limb)",
      label  = gsub("EtOH_","EtOH ",gsub("_"," ",colnames(etoh_6d))),
      stringsAsFactors = FALSE
    )
    etoh_df_p[[xcol]] <- etoh_6d[i, ]
    etoh_df_p[[ycol]] <- etoh_6d[j, ]

    # Nudge
    x_range <- diff(range(c(mds_df[[xcol]], bms_df_p[[xcol]], etoh_df_p[[xcol]])))
    nudge    <- 0.028 * x_range
    bms_true_x  <- bms_df_p[[xcol]]; bms_true_y  <- bms_df_p[[ycol]]
    etoh_true_x <- etoh_df_p[[xcol]]; etoh_true_y <- etoh_df_p[[ycol]]
    bms_df_p[[xcol]]  <- bms_df_p[[xcol]]  + nudge
    bms_df_p[[ycol]]  <- bms_df_p[[ycol]]  + nudge
    etoh_df_p[[xcol]] <- etoh_df_p[[xcol]] - nudge
    etoh_df_p[[ycol]] <- etoh_df_p[[ycol]] - nudge

    # Connector segments
    seg_df <- do.call(rbind, lapply(seq_len(nrow(cond_pairs)), function(k) {
      bi <- which(colnames(bms_6d)  == cond_pairs$bms_samp[k])
      ei <- which(colnames(etoh_6d) == cond_pairs$etoh_samp[k])
      if (length(bi) == 0 || length(ei) == 0) return(NULL)
      data.frame(x1 = bms_true_x[bi], y1 = bms_true_y[bi],
                 x2 = etoh_true_x[ei], y2 = etoh_true_y[ei])
    }))

    # Ring overlay data (reference coords only)
    ring_df <- mds_df[mds_df$sample %in% baseline_samps, ]

    ref_cols_needed <- c("sample","label","group", xcol, ycol)
    bms_cols_needed  <- c("sample","label","group", xcol, ycol)
    all_df <- rbind(mds_df[, ref_cols_needed], bms_df_p[, bms_cols_needed], etoh_df_p[, bms_cols_needed])

    p_bms <- make_scatter_bms(
      all_df, seg_df, ring_df[, c("sample","label","group", xcol, ycol)],
      xcol, ycol, pct_mds[i], pct_mds[j],
      title_str    = sprintf("Centered MDS — %s vs %s | %s (no tadpole) + BMS/EtOH", xcol, ycol, norm_label),
      subtitle_str = sprintf("Spearman centroid projection | %d common genes | %s | Xenopus tadpole excluded",
                             length(common), norm_label)
    )
    fname_bms <- file.path(FIG_DIR, sprintf("bms_mds_%s_PC%dvPC%d.png", suffix, i, j))
    ggsave(fname_bms, p_bms, width = 10, height = 8, dpi = 200, bg = "white")
    message("  Saved: ", basename(fname_bms))
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
                              title_str    = sprintf("Scree — Direct PCA | %s (no tadpole)", norm_label),
                              subtitle_str = sprintf("%d common genes × %d samples (Xenopus tadpole excluded)",
                                                     length(common), ncol(ref_n)))
  ggsave(file.path(FIG_DIR, paste0("pca_scree_", suffix, ".png")),
         p_scree_pca, width = 8, height = 5, dpi = 200, bg = "white")

  for (pair in pc_pairs) {
    i <- pair[1]; j <- pair[2]
    xcol <- paste0("PC", i); ycol <- paste0("PC", j)
    p <- make_scatter(
      pca_df, xcol, ycol, pct_pca[i], pct_pca[j],
      title_str    = sprintf("Direct PCA — %s vs %s | %s (no tadpole)", xcol, ycol, norm_label),
      subtitle_str = sprintf("PCA on pseudobulk | %d common genes | %s | Xenopus tadpole excluded",
                             length(common), norm_label)
    )
    fname <- file.path(FIG_DIR, sprintf("pca_%s_PC%dvPC%d.png", suffix, i, j))
    ggsave(fname, p, width = 9, height = 7, dpi = 200, bg = "white")
    message("  Saved: ", basename(fname))
  }
}

message("\nDone. ", length(list.files(FIG_DIR, "*.png")),
        " figures saved to ", FIG_DIR, "/")
