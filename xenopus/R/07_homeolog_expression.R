# =============================================================================
# 07 — Exercise 06: Homeolog expression (.L vs .S) in X. laevis
# Input:  results/xen_merged.rds
# Output: PDF figures in results/homeolog_exercise/
#
# Addresses all 5 open questions from nxr-2026/03-annotation/exercises/
#   06-homeolog-expression/README.md
# Featured pair: msx1 (limb patterning / blastema marker)
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

results_dir <- file.path("results", "homeolog_exercise")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

xen_merged <- readRDS(file.path("results", "xen_merged.rds"))
xen_merged <- JoinLayers(xen_merged)

all_genes <- rownames(xen_merged)
message("Genes in reference: ", length(all_genes),
        "  |  Cells: ", ncol(xen_merged))

# =============================================================================
# Q3 — CATALOGUE .L / .S / BOTH (establishes what we're working with)
# =============================================================================

message("\n--- Q3: Cataloguing homeolog pairs ---")

l_genes <- grep("\\.L(\\.[0-9]+)?$", all_genes, value = TRUE)
s_genes <- grep("\\.S(\\.[0-9]+)?$", all_genes, value = TRUE)

base_L <- sub("\\.L(\\.[0-9]+)?$", "", l_genes)
base_S <- sub("\\.S(\\.[0-9]+)?$", "", s_genes)

bases_both   <- intersect(unique(base_L), unique(base_S))
bases_L_only <- setdiff(unique(base_L), unique(base_S))
bases_S_only <- setdiff(unique(base_S), unique(base_L))
n_no_suffix  <- length(all_genes) - length(l_genes) - length(s_genes)

message("Both .L & .S: ", length(bases_both))
message(".L only:      ", length(bases_L_only))
message(".S only:      ", length(bases_S_only))
message("No suffix:    ", n_no_suffix)

p_catalog <- ggplot(
    data.frame(
      category = factor(c("Both .L & .S",".L only",".S only","No suffix"),
                        levels = c("Both .L & .S",".L only",".S only","No suffix")),
      n = c(length(bases_both), length(bases_L_only),
            length(bases_S_only), n_no_suffix)
    ),
    aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.65, alpha = 0.9, color = "white") +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Both .L & .S" = "#1565C0", ".L only" = "#2E7D32",
                                ".S only" = "#B71C1C", "No suffix" = "grey60")) +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.12))) +
  labs(title    = "Homeolog status — X. laevis reference",
       subtitle = paste0(scales::comma(length(all_genes)), " total genes"),
       x = NULL, y = "Number of genes") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

pdf(file.path(results_dir, "Q3_homeolog_catalog.pdf"), width = 7, height = 4.5)
print(p_catalog)
dev.off()

# =============================================================================
# BUILD CANONICAL PAIR TABLE
# =============================================================================

message("\nBuilding canonical pair table...")
mat <- GetAssayData(xen_merged, layer = "data")

pair_df <- bind_rows(lapply(bases_both, function(b) {
  l <- grep(paste0("^", b, "\\.L(\\.[0-9]+)?$"), all_genes, value = TRUE)
  s <- grep(paste0("^", b, "\\.S(\\.[0-9]+)?$"), all_genes, value = TRUE)
  if (length(l) > 1) l <- l[which.max(Matrix::rowSums(mat[l, , drop = FALSE]))]
  if (length(s) > 1) s <- s[which.max(Matrix::rowSums(mat[s, , drop = FALSE]))]
  data.frame(base = b, gene_L = l, gene_S = s, stringsAsFactors = FALSE)
}))

pair_df$mean_L   <- Matrix::rowMeans(mat[pair_df$gene_L, ])
pair_df$mean_S   <- Matrix::rowMeans(mat[pair_df$gene_S, ])
pair_df$var_mean <- apply(pair_df, 1, function(r) {
  var((as.numeric(mat[r["gene_L"], ]) + as.numeric(mat[r["gene_S"], ])) / 2)
})

top_pairs <- pair_df |> arrange(desc(var_mean)) |> head(2000)
message("Total pairs: ", nrow(pair_df), "  |  Top 2000 for correlation")

# =============================================================================
# Q2 — GLOBAL PEARSON CORRELATION DISTRIBUTION
# =============================================================================

message("\n--- Q2: Pearson correlations (top 2000 pairs) ---")

top_pairs$pearson_r <- sapply(seq_len(nrow(top_pairs)), function(i) {
  cor(as.numeric(mat[top_pairs$gene_L[i], ]),
      as.numeric(mat[top_pairs$gene_S[i], ]),
      method = "pearson")
})

frac_high <- mean(top_pairs$pearson_r > 0.7, na.rm = TRUE)
frac_low  <- mean(top_pairs$pearson_r < 0.3, na.rm = TRUE)
message(sprintf("r > 0.7: %.1f%%  |  r < 0.3: %.1f%%",
                frac_high * 100, frac_low * 100))

p_corr <- ggplot(top_pairs, aes(x = pearson_r)) +
  geom_histogram(bins = 60, fill = "#1565C0", alpha = 0.8, color = "white") +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "#B71C1C", linewidth = 0.8) +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "#2E7D32", linewidth = 0.8) +
  annotate("text", x = 0.15, y = Inf,
           label = sprintf("%.0f%% r < 0.3", frac_low * 100),
           vjust = 1.6, size = 3.5, color = "#B71C1C") +
  annotate("text", x = 0.85, y = Inf,
           label = sprintf("%.0f%% r > 0.7", frac_high * 100),
           vjust = 1.6, size = 3.5, color = "#2E7D32") +
  labs(title    = "Pearson r: .L vs .S expression across cells",
       subtitle = "Top 2,000 highest-variance homeolog pairs | X. laevis limb regeneration",
       x = "Pearson r (.L vs .S)", y = "Number of pairs") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

pdf(file.path(results_dir, "Q2_correlation_distribution.pdf"), width = 7, height = 4.5)
print(p_corr)
dev.off()

# =============================================================================
# Q1 — SINGLE PAIR SCATTER: msx1.L vs msx1.S
# msx1 chosen as a biologically relevant patterning gene in limb regeneration
# =============================================================================

message("\n--- Q1: Single pair scatter (msx1) ---")

msx1_row <- pair_df[tolower(pair_df$base) == "msx1", ]
if (nrow(msx1_row) == 0) {
  message("  msx1 not found as a pair — using highest expressed pair")
  msx1_row <- pair_df |> arrange(desc(mean_L + mean_S)) |> head(1)
}

gene_l <- msx1_row$gene_L[1]
gene_s <- msx1_row$gene_S[1]
r_val  <- cor(as.numeric(mat[gene_l, ]), as.numeric(mat[gene_s, ]))
message("  ", gene_l, " vs ", gene_s, "  |  r = ", round(r_val, 3))

set.seed(42)
n_plot <- min(5000, ncol(xen_merged))
idx    <- sample(ncol(xen_merged), n_plot)

scatter_df <- data.frame(
  L         = as.numeric(mat[gene_l, idx]),
  S         = as.numeric(mat[gene_s, idx]),
  cell_type = as.character(xen_merged$cell_type)[idx]
)

p_scatter <- ggplot(scatter_df, aes(x = L, y = S, color = cell_type)) +
  geom_point(size = 0.7, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.5) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Pearson r = %.3f", r_val),
           hjust = 1.1, vjust = -0.6, size = 3.5, color = "grey20") +
  labs(title    = paste0(gene_l, "  vs  ", gene_s),
       subtitle = paste0("n = ", scales::comma(n_plot),
                         " cells (subsampled); dashed line = y = x"),
       x = paste("Log-norm. expr.", gene_l),
       y = paste("Log-norm. expr.", gene_s),
       color = "Cell type") +
  theme_classic(base_size = 10) +
  theme(plot.title      = element_text(face = "bold"),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"))

pdf(file.path(results_dir, "Q1_msx1_scatter.pdf"), width = 8, height = 5.5)
print(p_scatter)
dev.off()

# =============================================================================
# Q4 — CELL-TYPE DIVERGENCE (subfunctionalization candidates)
# =============================================================================

message("\n--- Q4: Cell-type divergence ---")

ct_vec    <- as.character(xen_merged$cell_type)
ct_levels <- sort(unique(ct_vec))

top_pairs$ct_divergence <- sapply(seq_len(nrow(top_pairs)), function(i) {
  l_expr <- as.numeric(mat[top_pairs$gene_L[i], ])
  s_expr <- as.numeric(mat[top_pairs$gene_S[i], ])
  diffs  <- vapply(ct_levels, function(ct) {
    idx <- which(ct_vec == ct)
    if (length(idx) < 5) return(NA_real_)
    mean(l_expr[idx]) - mean(s_expr[idx])
  }, numeric(1))
  var(diffs, na.rm = TRUE)
})

top_div <- top_pairs |> arrange(desc(ct_divergence)) |> head(10)
message("Top 10 most cell-type-divergent pairs:")
print(top_div[, c("base","gene_L","gene_S","pearson_r","ct_divergence")])

heat_df <- bind_rows(lapply(seq_len(nrow(top_div)), function(i) {
  l_expr <- as.numeric(mat[top_div$gene_L[i], ])
  s_expr <- as.numeric(mat[top_div$gene_S[i], ])
  bind_rows(lapply(ct_levels, function(ct) {
    idx <- which(ct_vec == ct)
    if (length(idx) < 5) return(NULL)
    data.frame(pair      = top_div$base[i], cell_type = ct,
               L_mean    = mean(l_expr[idx]), S_mean = mean(s_expr[idx]),
               diff_LS   = mean(l_expr[idx]) - mean(s_expr[idx]))
  }))
}))

p_heat <- ggplot(heat_df, aes(x = cell_type, y = pair, fill = diff_LS)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#B71C1C", mid = "white", high = "#1565C0",
                       midpoint = 0, name = "Mean\n.L − .S") +
  labs(title    = "Cell-type divergence: mean .L − .S expression",
       subtitle = "Top 10 most divergent homeolog pairs — X. laevis",
       x = NULL, y = "Gene pair") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        plot.title  = element_text(face = "bold"))

pdf(file.path(results_dir, "Q4_divergence_heatmap.pdf"), width = 10, height = 5)
print(p_heat)
dev.off()

# Focused bar chart for top subfunctionalization candidate
top1     <- top_div[1, ]
focus_df <- heat_df |>
  filter(pair == top1$base) |>
  pivot_longer(c(L_mean, S_mean), names_to = "homeolog", values_to = "expr") |>
  mutate(homeolog = ifelse(homeolog == "L_mean", top1$gene_L, top1$gene_S))

p_focus <- ggplot(focus_df, aes(x = cell_type, y = expr, fill = homeolog)) +
  geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.9,
           color = "white") +
  scale_fill_manual(values = setNames(c("#1565C0","#B71C1C"),
                                       c(top1$gene_L, top1$gene_S))) +
  labs(title    = paste0("Subfunctionalization candidate: ", top1$base),
       subtitle = paste0(top1$gene_L, " (blue)  vs  ", top1$gene_S, " (red)"),
       x = NULL, y = "Mean log-norm. expr.", fill = "Homeolog") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
        plot.title  = element_text(face = "bold"))

pdf(file.path(results_dir, "Q4_top_subfunctionalized_pair.pdf"), width = 10, height = 4.5)
print(p_focus)
dev.off()

# =============================================================================
# SUMMARY
# =============================================================================

message("\n============================================================")
message("EXERCISE 06 — HOMEOLOG ANALYSIS SUMMARY")
message("------------------------------------------------------------")
message("Q3 | Both .L & .S: ", length(bases_both),
        "  | .L only: ", length(bases_L_only),
        "  | .S only: ", length(bases_S_only))
message("Q2 | r > 0.7 (co-expressed):     ", sprintf("%.1f%%", frac_high * 100))
message("   | r < 0.3 (divergent):        ", sprintf("%.1f%%", frac_low  * 100))
message("Q1 | msx1 Pearson r:             ", round(r_val, 3))
message("Q4 | Top subfunctionalization candidate: ", top1$base)
message("------------------------------------------------------------")
message("Q5 | RECOMMENDATION:")
message("   | If r > 0.7 dominates -> ambient RNA may inflate correlations.")
message("   | GOI with subfunctionalization evidence (e.g. msx1, sall1)")
message("   | should be kept as .L and .S separately for cell-type analyses.")
message("   | For bulk-like summaries, collapsing (.L + .S sum) is appropriate.")
message("============================================================")
message("Figures: ", results_dir)
