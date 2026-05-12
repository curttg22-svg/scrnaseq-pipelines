# =============================================================================
# plot_reamputation_spearman.R
#
# Single-row heatmap of Spearman rho for the Reamputation bulk RNA-seq
# condition vs all scRNA-seq pseudobulk timepoints across all three species.
#
# Reads from pre-computed CSVs — no RDS files needed.
# Run from comparison/ directory: source("R/plot_reamputation_spearman.R")
#
# Output: results/reamputation_spearman_heatmap.pdf
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

RESULTS_DIR <- "results"

# =============================================================================
# 1. LOAD AND RESHAPE
# =============================================================================

cor_mat <- read.csv(file.path(RESULTS_DIR, "bulk_scrna_correlation_matrix.csv"),
                    row.names = 1, check.names = FALSE)

df <- cor_mat["Reamputation", , drop = FALSE] |>
  as.data.frame() |>
  tibble::rownames_to_column("bulk_condition") |>
  pivot_longer(-bulk_condition, names_to = "timepoint", values_to = "rho")

# =============================================================================
# 2. ANNOTATE
# =============================================================================

df <- df |>
  mutate(
    prefix = sub("_.*", "", timepoint),
    species = case_when(
      prefix == "Leigh"   ~ "Axolotl\n(Leigh et al. 2018)",
      prefix == "LinBL"   ~ "Xenopus froglet\n(Lin et al. 2021)",
      prefix == "LinLBst" ~ "Xenopus tadpole\n(Lin et al. 2021)",
      prefix == "Mouse"   ~ "Mouse digit tip\n(Johnson, Masias & Lehoczky 2020)",
      TRUE ~ prefix
    ),
    tp_label = sub("^[^_]+_", "", timepoint)
  )

tp_order <- c(
  "Leigh_Intact", "Leigh_WH", "Leigh_EB", "Leigh_MB",
  "LinBL_0dpa", "LinBL_3dpa", "LinBL_7-14dpa", "LinBL_14dpa", "LinBL_14-52dpa",
  "LinLBst_NF50", "LinLBst_NF51", "LinLBst_NF52",
  "Mouse_0dpa", "Mouse_11dpa", "Mouse_12dpa", "Mouse_14dpa", "Mouse_17dpa"
)

df$timepoint <- factor(df$timepoint, levels = tp_order)
df$species   <- factor(df$species,
                       levels = c("Axolotl\n(Leigh et al. 2018)",
                                  "Xenopus froglet\n(Lin et al. 2021)",
                                  "Xenopus tadpole\n(Lin et al. 2021)",
                                  "Mouse digit tip\n(Johnson, Masias & Lehoczky 2020)"))

# Keep only the best-matching timepoint per species
df_best <- df |>
  group_by(species) |>
  slice_max(rho, n = 1, with_ties = FALSE) |>
  ungroup()

best_rho  <- max(df_best$rho)
rho_range <- range(df_best$rho)

# =============================================================================
# 3. PLOT
# =============================================================================

p <- ggplot(df_best, aes(x = timepoint, y = "Reamputation", fill = rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", rho),
                color = rho < mean(rho_range)),
            size = 3.2, fontface = "bold") +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#B2182B",
    midpoint = mean(rho_range),
    limits   = c(floor(min(df$rho) * 20) / 20,
                 ceiling(max(df$rho) * 20) / 20),
    name     = "Spearman rho"
  ) +
  scale_color_manual(values = c("TRUE" = "grey20", "FALSE" = "white"),
                     guide = "none") +
  scale_x_discrete(labels = function(x) sub("^[^_]+_", "", x)) +
  facet_grid(. ~ species, scales = "free_x", space = "free_x") +
  labs(
    title    = "Reamputation: bulk RNA-seq vs scRNA-seq pseudobulk (Spearman rho)",
    subtitle = sprintf(
      "Best match: Xenopus froglet 7-14 dpa (rho = %.3f)  |  Genes used: Axolotl 5,617  |  Xenopus 9,895  |  Mouse 10,987",
      best_rho
    ),
    x = "scRNA-seq timepoint",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 10),
    # strip prefix from x-axis tick labels (Leigh_WH -> WH, LinBL_0dpa -> 0dpa)
    axis.text.y      = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey93", color = NA),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 8.5, color = "grey40"),
    legend.position  = "right",
    legend.key.height = unit(1.2, "cm")
  )

out_file <- file.path(RESULTS_DIR, "reamputation_spearman_bestmatch.pdf")
pdf(out_file, width = 7, height = 3)
print(p)
dev.off()
message("Saved ", out_file)
