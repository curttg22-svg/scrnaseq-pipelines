# =============================================================================
# plot_regeneration_pca.R
#
# Builds a regenerative-spectrum reference from scRNA-seq pseudobulk
# (axolotl collab, Xenopus, mouse digit tip) and projects BMS/EtOH bulk
# axolotl limb RNA-seq onto it.
#
# Three normalizations are compared:
#   raw      — absolute pseudobulk expression (no normalization)
#   centered — subtract per-gene species mean before MDS
#   zscored  — subtract per-gene species mean and divide by species SD
#
# Output files (figures/ and results/) carry the suffix _centered / _zscored
# so all three can be compared side-by-side.
#
# Run from comparison/ directory: source("R/plot_regeneration_pca.R")
# =============================================================================

library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)

CACHE <- "results/pseudobulk_cache.rds"
dir.create("figures", showWarnings = FALSE)

# =============================================================================
# 1. Load cache — runs once
# =============================================================================

cache    <- readRDS(CACHE)
common   <- cache$common_genes
pseudos  <- cache$pseudos
bulk_mat <- cache$bulk_mat

axo_full     <- pseudos[["Axolotl (collaborator)"]]
mus_full     <- pseudos[["Mouse digit tip"]]
xen_all_full <- pseudos[["Xenopus (Lin 2021)"]]

axo_mat  <- axo_full[common, ]
mus_mat  <- mus_full[common, ]

xen_keep <- !colnames(xen_all_full) %in% c("LinBL_7-14dpa", "LinBL_14-52dpa")
xen_mat  <- xen_all_full[common, xen_keep, drop = FALSE]
xen_full <- xen_all_full[, xen_keep, drop = FALSE]
message("Xenopus: keeping ", sum(xen_keep), " of ", ncol(xen_all_full),
        " pseudobulk samples (excluded 2 pooled timepoints)")

# Raw common-gene reference matrix (9,254 genes × 15 samples)
ref_mat <- cbind(axo_mat, xen_mat, mus_mat)
message("Reference: ", ncol(ref_mat), " pseudobulk samples x ",
        nrow(ref_mat), " common genes")

# Full-resolution matrices for pairwise projection
ref_species_mats <- list(Axolotl = axo_full, Xenopus = xen_full, Mouse = mus_full)
ref_sample_species <- c(
  setNames(rep("Axolotl", ncol(axo_mat)), colnames(axo_mat)),
  setNames(rep("Xenopus", ncol(xen_mat)), colnames(xen_mat)),
  setNames(rep("Mouse",   ncol(mus_mat)), colnames(mus_mat))
)

bms_cols  <- grep("^BMS_|^Pretreatment$|^Reamputation$", colnames(bulk_mat), value = TRUE)
etoh_cols <- grep("^EtOH_", colnames(bulk_mat), value = TRUE)

# =============================================================================
# 2. Helper functions — defined once
# =============================================================================

# --- Normalization -----------------------------------------------------------
# Apply within-species centering or z-scoring to the common-gene ref matrix.
# Each gene's mean (center) or mean+SD (zscore) is computed across all samples
# of that species and subtracted / divided out, removing inter-species baseline
# differences so the MDS captures injury/regeneration dynamics rather than
# species identity.
normalize_ref_mat <- function(mat, species_vec, method) {
  if (method == "none") return(mat)
  result <- mat
  for (sp in unique(species_vec)) {
    cols   <- names(species_vec)[species_vec == sp]
    sp_mat <- mat[, cols, drop = FALSE]
    gm     <- rowMeans(sp_mat)
    if (method == "center") {
      result[, cols] <- sp_mat - gm
    } else {   # zscore
      gs <- apply(sp_mat, 1, sd); gs[gs < 1e-9] <- 1
      result[, cols] <- (sp_mat - gm) / gs
    }
  }
  result
}

# Apply the same normalization to the full-resolution species matrices used
# for pairwise projection. Each matrix contains one species only, so the
# within-species mean/SD is computed across all columns of that matrix.
normalize_species_mats <- function(mats_list, method) {
  if (method == "none") return(mats_list)
  lapply(mats_list, function(sp_mat) {
    gm <- rowMeans(sp_mat)
    if (method == "center") return(sp_mat - gm)
    gs <- apply(sp_mat, 1, sd); gs[gs < 1e-9] <- 1
    (sp_mat - gm) / gs
  })
}

# --- Reference metadata ------------------------------------------------------
make_meta <- function(coords) {
  nms <- rownames(coords)
  group <- dplyr::case_when(
    grepl("NonRegen", nms) & grepl("Axo",   nms) ~ "Axolotl (non-regen)",
    grepl("Regen",    nms) & grepl("Axo",   nms) ~ "Axolotl (regenerating)",
    grepl("^LinBL_",  nms)                        ~ "Xenopus froglet (non-regen)",
    grepl("^LinLBst_",nms)                        ~ "Xenopus tadpole (regen)",
    grepl("^Mouse_",  nms)                        ~ "Mouse digit tip",
    TRUE                                           ~ "Unknown"
  )
  label <- nms |>
    gsub("AxoCollab_", "",      x = _) |>
    gsub("LinBL_",     "BL ",   x = _) |>
    gsub("LinLBst_",   "LBst ", x = _) |>
    gsub("Mouse_",     "Mouse ",x = _) |>
    gsub("_",          " ",     x = _)
  data.frame(sample = nms, label = label, group = group,
             MDS1 = coords[, 1], MDS2 = coords[, 2],
             stringsAsFactors = FALSE)
}

# --- Pairwise projection onto diagonal axis ----------------------------------
# For each reference sample r (species S), rho is computed over the genes
# present in both S's full-resolution matrix and the new sample — maximising
# gene count per comparison without requiring a 3-way common set.
project_to_linear <- function(new_mat, ref_species_mats, ref_sample_species,
                               ref_coords, linear_pca, score_ranked, group_label) {
  ref_samples <- rownames(ref_coords)
  n_new       <- ncol(new_mat)
  rho <- matrix(0, nrow = length(ref_samples), ncol = n_new,
                dimnames = list(ref_samples, colnames(new_mat)))
  for (r in ref_samples) {
    sp_mat     <- ref_species_mats[[ref_sample_species[r]]]
    pair_genes <- intersect(rownames(sp_mat), rownames(new_mat))
    ref_vec    <- sp_mat[pair_genes, r, drop = FALSE]
    new_sub    <- new_mat[pair_genes, , drop = FALSE]
    rho[r, ]   <- cor(ref_vec, new_sub, method = "spearman")[1, ]
  }
  centroid_2d <- vapply(colnames(new_mat), function(s) {
    w <- pmax(rho[, s], 0)
    if (sum(w) == 0) return(c(0, 0))
    colSums(ref_coords * w) / sum(w)
  }, numeric(2))
  centered    <- t(centroid_2d) - matrix(linear_pca$center, nrow = n_new, ncol = 2, byrow = TRUE)
  proj_scores <- centered %*% linear_pca$rotation
  raw_score   <- proj_scores[, 1]
  raw_perp    <- proj_scores[, 2]
  ref_raw     <- sort(linear_pca$x[, 1])
  ref_ranks   <- score_ranked[order(linear_pca$x[, 1])]
  ranked_new  <- approx(ref_raw, ref_ranks, xout = raw_score, rule = 2)$y
  perp_new    <- (raw_perp / max(abs(linear_pca$x[, 2]))) * 0.30
  label <- colnames(new_mat) |>
    gsub("BMS_",  "BMS ",  x = _) |>
    gsub("EtOH_", "EtOH ", x = _) |>
    gsub("_",     " ",     x = _)
  data.frame(
    sample  = colnames(new_mat), label = label, group = group_label,
    score   = ranked_new,
    spine_x = 1 - ranked_new, spine_y = ranked_new,
    x_plot  = (1 - ranked_new) + perp_new / sqrt(2),
    y_plot  = ranked_new       + perp_new / sqrt(2),
    stringsAsFactors = FALSE
  )
}

# --- Plot style constants ----------------------------------------------------
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
all_colors <- c(group_colors,
                "BMS (axolotl limb)"  = "#FDD835",
                "EtOH (axolotl limb)" = "#E53935")
all_shapes <- c(group_shapes,
                "BMS (axolotl limb)"  = 23,
                "EtOH (axolotl limb)" = 23)

make_diagonal_plot <- function(diag_df, subtitle_text) {
  pad      <- 0.06
  s_seq    <- seq(0 - pad, 1 + pad, length.out = 150)
  backbone <- data.frame(x = 1 - s_seq, y = s_seq)
  grad_df  <- data.frame(x = 1 - s_seq, y = s_seq, fill = s_seq)
  ggplot(diag_df, aes(x_plot, y_plot)) +
    geom_tile(data = grad_df, aes(x = x, y = y, fill = fill),
              width = 0.55, height = 0.55, alpha = 0.40, inherit.aes = FALSE) +
    scale_fill_gradient(low = "#7B0000", high = "#0D2B6B", guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_line(data = backbone, aes(x = x, y = y),
              color = "grey55", linewidth = 1.0, inherit.aes = FALSE) +
    geom_segment(aes(xend = spine_x, yend = spine_y),
                 color = "grey55", linewidth = 0.3, alpha = 0.7) +
    geom_point(aes(fill = group, shape = group),
               size = 6, color = "white", stroke = 0.5, alpha = 0.97) +
    geom_text_repel(aes(label = label), size = 3.0, color = "grey92",
                    bg.color = "#0d0d1a", bg.r = 0.12,
                    box.padding = 0.5, point.padding = 0.3,
                    max.overlaps = 30,
                    segment.color = "grey55", segment.size = 0.3) +
    annotate("text", x = 1 + pad, y = 0 - pad,
             label = "Less regenerative", color = "#EF9A9A",
             hjust = 1, vjust = 1.8, size = 3.5, fontface = "italic") +
    annotate("text", x = 0 - pad, y = 1 + pad,
             label = "More regenerative", color = "#90CAF9",
             hjust = 0, vjust = -0.8, size = 3.5, fontface = "italic") +
    scale_fill_manual(values = all_colors,  name = NULL) +
    scale_shape_manual(values = all_shapes, name = NULL) +
    coord_fixed(ratio = 1,
                xlim = range(diag_df$x_plot) + c(-1, 1) * 0.12,
                ylim = range(diag_df$y_plot) + c(-1, 1) * 0.12) +
    labs(title = "Regenerative spectrum", subtitle = subtitle_text,
         x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      plot.background   = element_rect(fill = "#0d0d1a", color = NA),
      panel.background  = element_rect(fill = "#0d0d1a", color = NA),
      panel.grid        = element_blank(),
      axis.text         = element_blank(), axis.ticks = element_blank(),
      plot.title        = element_text(color = "grey97", face = "bold", size = 16,
                                       margin = margin(b = 4)),
      plot.subtitle     = element_text(color = "grey55", size = 9),
      legend.text       = element_text(color = "grey85", size = 10),
      legend.key        = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "#0d0d1a", color = NA),
      legend.position   = "right",
      plot.margin       = margin(20, 20, 20, 20)
    ) +
    guides(fill  = guide_legend(override.aes = list(size = 5)),
           shape = guide_legend(override.aes = list(size = 5)))
}

make_inset_plot <- function(bulk_all_df) {
  n_bi        <- nrow(bulk_all_df)
  within_rank <- rank(bulk_all_df$score, ties.method = "average")
  inset_score <- (within_rank - 1) / (n_bi - 1)
  perp_raw    <- bulk_all_df$x_plot - bulk_all_df$spine_x
  max_pr      <- max(abs(perp_raw))
  perp_i      <- if (max_pr > 1e-9) perp_raw / max_pr * 0.25 else rep(0, n_bi)
  df <- data.frame(
    sample  = bulk_all_df$sample, label = bulk_all_df$label, group = bulk_all_df$group,
    score   = inset_score,
    spine_x = 1 - inset_score, spine_y = inset_score,
    x_plot  = (1 - inset_score) + perp_i / sqrt(2),
    y_plot  = inset_score       + perp_i / sqrt(2)
  )
  pad_i  <- 0.06
  s_i    <- seq(0 - pad_i, 1 + pad_i, length.out = 120)
  bb_i   <- data.frame(x = 1 - s_i, y = s_i)
  gr_i   <- data.frame(x = 1 - s_i, y = s_i, fill = s_i)
  spread <- diff(range(bulk_all_df$score))
  ggplot(df, aes(x_plot, y_plot)) +
    geom_tile(data = gr_i, aes(x = x, y = y, fill = fill),
              width = 0.55, height = 0.55, alpha = 0.35, inherit.aes = FALSE) +
    scale_fill_gradient(low = "#7B0000", high = "#0D2B6B", guide = "none") +
    ggnewscale::new_scale_fill() +
    geom_line(data = bb_i, aes(x = x, y = y),
              color = "grey55", linewidth = 0.8, inherit.aes = FALSE) +
    geom_segment(aes(xend = spine_x, yend = spine_y),
                 color = "grey55", linewidth = 0.25, alpha = 0.7) +
    geom_point(aes(fill = group, shape = group),
               size = 5, color = "white", stroke = 0.5, alpha = 0.97) +
    geom_text_repel(aes(label = label), size = 2.8, color = "grey92",
                    bg.color = "#0d0d1a", bg.r = 0.10,
                    box.padding = 0.45, point.padding = 0.25,
                    max.overlaps = 20,
                    segment.color = "grey55", segment.size = 0.25) +
    annotate("text", x = 1 + pad_i, y = 0 - pad_i,
             label = "Less\nregenerative", color = "#EF9A9A",
             hjust = 1, vjust = 1.5, size = 2.4, fontface = "italic") +
    annotate("text", x = 0 - pad_i, y = 1 + pad_i,
             label = "More\nregenerative", color = "#90CAF9",
             hjust = 0, vjust = -0.5, size = 2.4, fontface = "italic") +
    scale_fill_manual(values  = all_colors[c("BMS (axolotl limb)", "EtOH (axolotl limb)")],
                      name = NULL) +
    scale_shape_manual(values = all_shapes[c("BMS (axolotl limb)", "EtOH (axolotl limb)")],
                       name = NULL) +
    coord_fixed(ratio = 1,
                xlim = range(df$x_plot) + c(-1, 1) * 0.15,
                ylim = range(df$y_plot) + c(-1, 1) * 0.15) +
    labs(title    = "BMS / EtOH conditions — magnified",
         subtitle = sprintf("Within-cluster relative ordering  |  actual score spread: %.4f", spread)) +
    theme_minimal(base_size = 11) +
    theme(
      plot.background   = element_rect(fill = "#0d0d1a", color = "grey50", linewidth = 0.8),
      panel.background  = element_rect(fill = "#0d0d1a", color = NA),
      panel.grid        = element_blank(),
      axis.text         = element_blank(), axis.ticks = element_blank(),
      plot.title        = element_text(color = "grey95", face = "bold", size = 11,
                                       margin = margin(b = 2)),
      plot.subtitle     = element_text(color = "grey55", size = 7.5),
      legend.text       = element_text(color = "grey85", size = 8),
      legend.key        = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "#0d0d1a", color = NA),
      legend.position   = "bottom",
      plot.margin       = margin(8, 8, 8, 8)
    ) +
    guides(fill  = guide_legend(override.aes = list(size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))
}

# =============================================================================
# 3. Main analysis function — called once per normalization method
# =============================================================================

run_spectrum_analysis <- function(ref_mat_n, ref_species_mats_n, suffix, norm_label) {

  message("\n========== Normalization: ", norm_label, " ==========")

  # --- Reference MDS ---------------------------------------------------------
  message("Computing Spearman MDS...")
  rho_ref    <- cor(ref_mat_n, method = "spearman")
  mds        <- cmdscale(as.dist(1 - rho_ref), k = 2, eig = TRUE)
  ref_coords <- mds$points
  eig_pos    <- pmax(mds$eig, 0)
  pct_var    <- round(100 * eig_pos[1:2] / sum(eig_pos), 1)
  message("MDS1: ", pct_var[1], "%   MDS2: ", pct_var[2], "%")

  ref_meta <- make_meta(ref_coords)

  # --- Reference map (2D MDS scatter) ----------------------------------------
  p_map <- ggplot(ref_meta, aes(MDS1, MDS2, fill = group, shape = group, label = label)) +
    stat_ellipse(aes(color = group, group = group),
                 type = "norm", level = 0.8,
                 linewidth = 0.4, linetype = "dashed", alpha = 0.5,
                 show.legend = FALSE) +
    geom_point(size = 4.5, color = "white", stroke = 0.4, alpha = 0.95) +
    geom_text_repel(size = 3, color = "grey90", bg.color = "grey10", bg.r = 0.12,
                    box.padding = 0.45, point.padding = 0.3,
                    max.overlaps = 25, segment.color = "grey50", segment.size = 0.3) +
    scale_fill_manual(values  = group_colors, name = NULL) +
    scale_color_manual(values = group_colors, name = NULL) +
    scale_shape_manual(values = group_shapes, name = NULL) +
    labs(
      title    = paste0("Regenerative spectrum — reference map (", norm_label, ")"),
      subtitle = paste0("Spearman MDS  |  ", length(common), " common genes  |  ",
                        "MDS1 ", pct_var[1], "%  MDS2 ", pct_var[2], "%"),
      x = "MDS1", y = "MDS2"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "#1a1a2e", color = NA),
      panel.background = element_rect(fill = "#1a1a2e", color = NA),
      panel.grid.major = element_line(color = "#2e2e4e", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(color = "grey95", face = "bold", size = 14),
      plot.subtitle    = element_text(color = "grey65", size = 9),
      axis.title       = element_text(color = "grey80"),
      axis.text        = element_text(color = "grey60"),
      legend.text      = element_text(color = "grey85", size = 9),
      legend.key       = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = "#1a1a2e", color = NA),
      plot.margin      = margin(16, 16, 16, 16)
    ) +
    guides(fill  = guide_legend(override.aes = list(size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))

  ggsave(paste0("results/regeneration_pca", suffix, ".pdf"),
         p_map, width = 10, height = 7.5)
  ggsave(paste0("figures/regeneration_pca", suffix, ".png"),
         p_map, width = 10, height = 7.5, dpi = 200, bg = "#1a1a2e")
  message("Saved: figures/regeneration_pca", suffix, ".png")

  # --- Diagonal linearization ------------------------------------------------
  linear_pca   <- prcomp(ref_coords, center = TRUE, scale. = FALSE)
  score        <- linear_pca$x[, 1]
  perp         <- linear_pca$x[, 2]
  n_samp       <- length(score)
  rank_pos     <- rank(score, ties.method = "average")
  score_ranked <- (rank_pos - 1) / (n_samp - 1)

  axo_idx <- ref_meta$group == "Axolotl (regenerating)"
  if (mean(rank_pos[axo_idx]) < n_samp / 2) score_ranked <- 1 - score_ranked

  perp_scaled <- (perp / max(abs(perp))) * 0.30

  diag_cols <- c("sample", "label", "group", "score",
                 "spine_x", "spine_y", "x_plot", "y_plot")

  ref_diag          <- ref_meta
  ref_diag$score    <- score_ranked
  ref_diag$spine_x  <- 1 - score_ranked
  ref_diag$spine_y  <- score_ranked
  ref_diag$x_plot   <- ref_diag$spine_x + perp_scaled / sqrt(2)
  ref_diag$y_plot   <- ref_diag$spine_y + perp_scaled / sqrt(2)

  # --- Bulk projection -------------------------------------------------------
  bulk_bms  <- project_to_linear(bulk_mat[, bms_cols,  drop = FALSE],
                                  ref_species_mats_n, ref_sample_species,
                                  ref_coords, linear_pca, score_ranked,
                                  group_label = "BMS (axolotl limb)")
  bulk_etoh <- project_to_linear(bulk_mat[, etoh_cols, drop = FALSE],
                                  ref_species_mats_n, ref_sample_species,
                                  ref_coords, linear_pca, score_ranked,
                                  group_label = "EtOH (axolotl limb)")

  bulk_all_df <- rbind(bulk_bms, bulk_etoh)
  diag_df     <- rbind(ref_diag[, diag_cols],
                       bulk_bms[,  diag_cols],
                       bulk_etoh[, diag_cols])

  message(sprintf("BMS/EtOH score range: %.4f – %.4f  (spread %.4f)",
                  min(bulk_all_df$score), max(bulk_all_df$score),
                  diff(range(bulk_all_df$score))))

  # --- Main diagonal plot ----------------------------------------------------
  pct_lp <- round(100 * linear_pca$sdev[1]^2 / sum(linear_pca$sdev^2), 1)
  sub_main <- paste0("PC1 of Spearman MDS  |  ", length(common),
                     " common genes  |  ", ncol(ref_mat_n), " reference samples  |  ",
                     norm_label, "  |  ", pct_lp, "% of MDS variance")

  cx <- mean(bulk_all_df$x_plot); cy <- mean(bulk_all_df$y_plot); box_r <- 0.045
  p_diag <- make_diagonal_plot(diag_df, sub_main) +
    annotate("rect",
             xmin = cx - box_r, xmax = cx + box_r,
             ymin = cy - box_r, ymax = cy + box_r,
             color = "white", fill = NA, linetype = "dashed",
             linewidth = 0.7, alpha = 0.8) +
    annotate("text", x = cx + box_r + 0.01, y = cy + box_r + 0.01,
             label = "BMS / EtOH\n(see inset)", color = "grey80",
             hjust = 0, vjust = 0, size = 2.8, fontface = "italic")

  ggsave(paste0("results/regeneration_linear", suffix, ".pdf"),
         p_diag, width = 13, height = 8)
  ggsave(paste0("figures/regeneration_linear", suffix, ".png"),
         p_diag, width = 13, height = 8, dpi = 200, bg = "#0d0d1a")
  message("Saved: figures/regeneration_linear", suffix, ".png")

  # --- Inset + combined ------------------------------------------------------
  p_inset <- make_inset_plot(bulk_all_df)
  ggsave(paste0("figures/regeneration_linear", suffix, "_inset.png"),
         p_inset, width = 6, height = 5, dpi = 200, bg = "#0d0d1a")
  message("Saved: figures/regeneration_linear", suffix, "_inset.png")

  p_combined <- p_diag +
    inset_element(p_inset, left = 0.00, bottom = 0.00,
                  right = 0.46, top = 0.50, align_to = "plot")
  ggsave(paste0("results/regeneration_linear", suffix, "_combined.pdf"),
         p_combined, width = 16, height = 9)
  ggsave(paste0("figures/regeneration_linear", suffix, "_combined.png"),
         p_combined, width = 16, height = 9, dpi = 200, bg = "#0d0d1a")
  message("Saved: figures/regeneration_linear", suffix, "_combined.png")

  invisible(list(ref_coords = ref_coords, linear_pca = linear_pca,
                 score_ranked = score_ranked, diag_df = diag_df))
}

# =============================================================================
# 4. Run all three normalizations
# =============================================================================

for (method in c("none", "center", "zscore")) {
  suffix     <- switch(method, none = "", center = "_centered", zscore = "_zscored")
  norm_label <- switch(method,
                       none   = "raw (no normalization)",
                       center = "within-species centered",
                       zscore = "within-species z-scored")

  ref_mat_n          <- normalize_ref_mat(ref_mat, ref_sample_species, method)
  ref_species_mats_n <- normalize_species_mats(ref_species_mats, method)

  run_spectrum_analysis(ref_mat_n, ref_species_mats_n, suffix, norm_label)
}

message("\nAll three normalizations complete.")

# =============================================================================
# 5. Regenerating-only reference — same three normalizations, regen samples only
#
# Excludes:
#   AxoCollab_NonRegen_limb  — intact non-regenerating axolotl limb
#   LinBL_*                  — Xenopus froglet (non-regenerative adult)
# Keeps:
#   AxoCollab_Regen_*        — axolotl blastema timepoints
#   LinLBst_*                — Xenopus tadpole (regenerative stage)
#   Mouse_*                  — mouse digit tip (regenerating)
# =============================================================================

regen_keep    <- !grepl("NonRegen|^LinBL_", colnames(ref_mat))
ref_mat_regen <- ref_mat[, regen_keep, drop = FALSE]
regen_sp_vec  <- ref_sample_species[colnames(ref_mat_regen)]

message("\nRegenerating-only reference: ", ncol(ref_mat_regen), " samples")
message("  ", paste(colnames(ref_mat_regen), collapse = "\n  "))

for (method in c("none", "center", "zscore")) {
  suffix_base <- switch(method, none = "", center = "_centered", zscore = "_zscored")
  suffix      <- paste0("_regen", suffix_base)
  norm_label  <- switch(method,
                        none   = "raw (no normalization)",
                        center = "within-species centered",
                        zscore = "within-species z-scored")

  ref_mat_n          <- normalize_ref_mat(ref_mat_regen, regen_sp_vec, method)
  ref_species_mats_n <- normalize_species_mats(ref_species_mats, method)

  run_spectrum_analysis(ref_mat_n, ref_species_mats_n, suffix, norm_label)
}

message("\nRegeneration-only analyses complete.")
message("  _zscored    — within-species z-scored")
