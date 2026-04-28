# =============================================================================
# Bulk RNA-seq vs scRNA-seq pseudobulk correlation — multi-species
#
# Run this script from the comparison/ directory (working directory matters
# for the source() call and results/ output path).
#
# Inputs:
#   Bulk (axolotl whole blastema, normalized counts):
#     ~/Desktop/All RNA normlized count csv files/BMS_normalized_counts.csv
#     ~/Desktop/All RNA normlized count csv files/24 hrs_normalized_counts.csv
#     ~/Desktop/All RNA normlized count csv files/Pretreatment_normalized_counts.csv
#     ~/Desktop/All RNA normlized count csv files/Reamputation_normalized_counts.csv
#
#   scRNA-seq datasets (one block per species below):
#     Axolotl (Leigh 2018): ~/Desktop/Axolotl_scrnaseq_project_V3/results/axo_integrated.rds
#     Xenopus (Lin 2021):   ~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds
#     Mouse digit tip:      set MOUSE_RDS and uncomment the mouse block (Section 3)
#
# Outputs:
#   results/bulk_scrna_correlation_matrix.csv
#   results/bulk_scrna_correlation_heatmap.pdf
#
# Method:
#   Each scRNA-seq dataset is collapsed to a genes x timepoints pseudobulk
#   matrix (mean log-normalized expression per timepoint). Gene names are
#   harmonized to HGNC uppercase symbols via utils_harmonize.R:
#     - Axolotl (Leigh 2018): org.Hs.eg.db ALIAS lookup (SwissProt -> HGNC)
#     - Xenopus: .L/.S suffix stripped, then toupper()
#     - Mouse: toupper() (mouse symbols are case-only variants of human)
#   Spearman correlation is computed between each bulk condition and each
#   scRNA-seq timepoint using genes shared across the two datasets.
#   LOC IDs are excluded from the bulk (uncharacterized genes with no
#   cross-dataset matches). Trinity IDs are excluded from Leigh 2018
#   (no HGNC mapping without sequence realignment to the genome).
#
# To add the mouse digit-tip dataset:
#   1. Set MOUSE_RDS to the path of the Seurat .rds file
#   2. Uncomment the mouse block in Section 3
#   3. Adjust tp_var / tp_order to match the object's metadata
#
# Dependencies:
#   Seurat, Matrix, ggplot2, dplyr, tidyr, AnnotationDbi
#   org.Hs.eg.db (for axolotl alias harmonization)
#   Install: BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db"))
# =============================================================================

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)

source("R/utils_harmonize.R")

# =============================================================================
# PATHS & CONFIG
# =============================================================================

BULK_DIR  <- "~/Desktop/All RNA normlized count csv files"
AXO_RDS   <- "~/Desktop/Axolotl_scrnaseq_project_V3/results/axo_integrated.rds"
XEN_RDS   <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
# MOUSE_RDS <- "~/path/to/mouse_digit.rds"   # set when dataset is available

RESULTS_DIR <- "results"
dir.create(RESULTS_DIR, showWarnings = FALSE)

WH_REPS     <- c("WH_N4", "WH_N5", "WH_N6")
BULK_ORDER  <- c("BMS_acute", "BMS_24hrs", "Pretreatment", "Reamputation")

# =============================================================================
# 1. BULK DATA
# =============================================================================

message("Loading bulk RNA-seq data...")

.load_bulk_condition <- function(path, condition) {
  df    <- read.csv(path, row.names = 1, check.names = FALSE)
  named <- df[!grepl("^LOC", rownames(df)), , drop = FALSE]
  setNames(data.frame(rowMeans(named)), condition)
}

bulk_list <- list(
  .load_bulk_condition(file.path(BULK_DIR, "BMS_normalized_counts.csv"),          "BMS_acute"),
  .load_bulk_condition(file.path(BULK_DIR, "24 hrs_normalized_counts.csv"),        "BMS_24hrs"),
  .load_bulk_condition(file.path(BULK_DIR, "Pretreatment_normalized_counts.csv"),  "Pretreatment"),
  .load_bulk_condition(file.path(BULK_DIR, "Reamputation_normalized_counts.csv"),  "Reamputation")
)

# Intersect genes across all four conditions, then cbind
bulk_genes <- Reduce(intersect, lapply(bulk_list, rownames))
bulk_raw   <- do.call(cbind, lapply(bulk_list, function(x) x[bulk_genes, , drop = FALSE]))

# Harmonize bulk gene names to canonical HGNC (alias pass catches edge cases
# where older symbols are used; official symbols map to themselves)
bulk_mat <- aggregate_to_hgnc(bulk_raw, species = "human")
message("Bulk: ", nrow(bulk_mat), " HGNC genes x ", ncol(bulk_mat), " conditions")

# =============================================================================
# 2. PSEUDOBULK HELPER
# =============================================================================

# Compute mean log-normalized expression per timepoint from a Seurat object.
# Returns a dense matrix: genes x timepoints (colnames prefixed by prefix_).
.mean_per_tp <- function(expr_mat, tp_factor, tp_order, prefix) {
  parts <- lapply(tp_order, function(tp) {
    idx <- which(tp_factor == tp)
    if (length(idx) == 0L) return(NULL)
    m <- Matrix::rowMeans(expr_mat[, idx, drop = FALSE])
    # data.frame(m) preserves names(m) as rownames; as.numeric() would strip them
    setNames(data.frame(m), paste0(prefix, "_", tp))
  })
  do.call(cbind, Filter(Negate(is.null), parts))
}

# =============================================================================
# 3. DATASET-SPECIFIC PSEUDOBULK FUNCTIONS
#
# Each function returns a dense matrix:
#   rows    = HGNC-harmonized gene symbols (deduplicated by averaging)
#   columns = prefixed timepoint labels (e.g., "Leigh_Intact", "Lin_0dpa")
# =============================================================================

# --- Axolotl Leigh 2018 (Trinity assembly, SwissProt-annotated names) -------

pseudobulk_axolotl_leigh <- function(rds_path) {
  message("Computing pseudobulk: Axolotl (Leigh 2018)...")
  axo    <- JoinLayers(readRDS(rds_path))
  expr   <- GetAssayData(axo, layer = "data")

  # Pool WH_N4/N5/N6 replicates into a single "WH" timepoint
  tp           <- as.character(axo$timepoint)
  tp[tp %in% WH_REPS] <- "WH"
  tp_order     <- c("Intact", "WH", "EB", "MB")
  tp_factor    <- factor(tp, levels = tp_order)

  pb <- .mean_per_tp(expr, tp_factor, tp_order, prefix = "Leigh")
  message("  Axolotl: ", ncol(expr), " cells x ", nrow(expr), " genes")

  # Exclude Trinity IDs (^c[0-9]) — no HGNC mapping without sequence alignment.
  # Harmonize remaining named genes via org.Hs.eg.db ALIAS table:
  #   catches NCBI-listed synonyms and some SwissProt abbreviations.
  result <- aggregate_to_hgnc(
    as.matrix(pb),
    species   = "human",
    filter_fn = function(g) !grepl("^c[0-9]", g)
  )
  message("  After harmonization: ", nrow(result), " HGNC symbols retained")
  result
}

# --- Xenopus laevis Lin 2021 (pseudotetraploid, .L/.S homeologs) ------------

pseudobulk_xenopus_lin <- function(rds_path) {
  message("Computing pseudobulk: Xenopus (Lin 2021)...")
  xen    <- JoinLayers(readRDS(rds_path))
  expr   <- GetAssayData(xen, layer = "data")
  tp_order <- c("0dpa", "3dpa", "7-14dpa", "14dpa", "14-52dpa")
  tp_factor <- xen$timepoint

  # Sum .L and .S homeologs into one base-gene symbol per cell using a sparse
  # indicator matrix multiply. Much faster than row-by-row aggregation.
  xen_genes    <- rownames(expr)
  base_symbols <- toupper(sub("\\.(L|S)(\\.[0-9]+)?$", "", xen_genes))
  unique_bases <- unique(base_symbols)
  message("  Xenopus: ", length(unique_bases), " base genes from ",
          length(xen_genes), " homeolog rows")

  indicator <- sparseMatrix(
    i = seq_along(xen_genes),
    j = match(base_symbols, unique_bases),
    x = 1,
    dims     = c(length(xen_genes), length(unique_bases)),
    dimnames = list(xen_genes, unique_bases)
  )
  expr_summed <- t(indicator) %*% expr   # base genes x cells

  pb <- .mean_per_tp(expr_summed, tp_factor, tp_order, prefix = "Lin")
  as.matrix(pb)   # rownames already HGNC-compatible uppercase base symbols
}

# --- Mouse digit tip (placeholder — fill in when dataset is available) ------
#
# pseudobulk_mouse_digit <- function(rds_path,
#                                     tp_var   = "timepoint",
#                                     tp_order = NULL,
#                                     prefix   = "Mouse") {
#   message("Computing pseudobulk: Mouse digit tip...")
#   mus      <- JoinLayers(readRDS(rds_path))
#   expr     <- GetAssayData(mus, layer = "data")
#   tp_vals  <- mus[[tp_var, drop = TRUE]]
#   if (is.null(tp_order)) tp_order <- levels(factor(tp_vals))
#   tp_factor <- factor(tp_vals, levels = tp_order)
#
#   pb <- .mean_per_tp(expr, tp_factor, tp_order, prefix = prefix)
#
#   # toupper() handles ~95%+ of mouse-human symbol pairs (same symbol,
#   # different capitalization). Genuine symbol divergences fall out at
#   # the intersection step — acceptable for transcriptome-wide Spearman.
#   aggregate_to_hgnc(as.matrix(pb), species = "mouse")
# }

# =============================================================================
# 4. BUILD PSEUDOBULK LIST
#    Each entry: HGNC-gene-symbol matrix x timepoints
# =============================================================================

pseudos <- list(
  "Axolotl (Leigh 2018)" = pseudobulk_axolotl_leigh(AXO_RDS),
  "Xenopus (Lin 2021)"   = pseudobulk_xenopus_lin(XEN_RDS)
  # Add mouse dataset when available:
  # "Mouse digit tip" = pseudobulk_mouse_digit(
  #   MOUSE_RDS,
  #   tp_var   = "timepoint",               # adjust to match metadata column
  #   tp_order = c("Intact","3dpa","7dpa","14dpa"),  # adjust to actual levels
  #   prefix   = "Mouse"
  # )
)

# =============================================================================
# 5. SPEARMAN CORRELATIONS
# =============================================================================

.spearman_cor_mat <- function(bulk_mat, scrna_mat) {
  shared <- intersect(rownames(bulk_mat), rownames(scrna_mat))
  b      <- bulk_mat[shared, , drop = FALSE]
  s      <- scrna_mat[shared, , drop = FALSE]
  cor_out <- matrix(
    NA_real_, nrow = ncol(b), ncol = ncol(s),
    dimnames = list(colnames(b), colnames(s))
  )
  for (i in colnames(b))
    for (j in colnames(s))
      cor_out[i, j] <- cor(b[, i], s[, j],
                           method = "spearman", use = "complete.obs")
  cor_out
}

message("\nComputing Spearman correlations...")
cor_list <- lapply(names(pseudos), function(nm) {
  shared_n <- length(intersect(rownames(bulk_mat), rownames(pseudos[[nm]])))
  message("  ", nm, ": ", shared_n, " shared genes")
  .spearman_cor_mat(bulk_mat, pseudos[[nm]])
})
names(cor_list) <- names(pseudos)

cor_all <- do.call(cbind, cor_list)

message("\nCorrelation matrix:")
print(round(cor_all, 3))

write.csv(round(cor_all, 4),
          file.path(RESULTS_DIR, "bulk_scrna_correlation_matrix.csv"))
message("Saved bulk_scrna_correlation_matrix.csv")

# =============================================================================
# 6. HEATMAP
# =============================================================================

# Derive species labels from column name prefixes (Axo_, Xen_, Mouse_, ...)
.prefix_to_species <- function(col_names) {
  prefixes <- sub("_.*", "", col_names)
  species_labels <- c(
    Leigh  = "Axolotl (Leigh 2018)",
    Gerber = "Axolotl (Gerber 2018)",
    Lin    = "Xenopus laevis (Lin 2021)",
    Mouse  = "Mouse digit tip"
  )
  ifelse(prefixes %in% names(species_labels),
         species_labels[prefixes], prefixes)
}

cor_df <- as.data.frame(cor_all) |>
  tibble::rownames_to_column("bulk_condition") |>
  pivot_longer(-bulk_condition,
               names_to  = "scrna_timepoint",
               values_to = "rho")

cor_df$bulk_condition  <- factor(cor_df$bulk_condition,  levels = BULK_ORDER)
cor_df$scrna_timepoint <- factor(cor_df$scrna_timepoint, levels = colnames(cor_all))
cor_df$species         <- .prefix_to_species(as.character(cor_df$scrna_timepoint))
cor_df$species         <- factor(cor_df$species,
                                  levels = unique(.prefix_to_species(colnames(cor_all))))

# Gene counts per dataset for subtitle
gene_counts <- sapply(names(pseudos), function(nm) {
  length(intersect(rownames(bulk_mat), rownames(pseudos[[nm]])))
})
subtitle_parts <- paste0(names(pseudos), ": ", gene_counts, " shared genes")

p <- ggplot(cor_df, aes(x = scrna_timepoint, y = bulk_condition, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 3, color = "black") +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#B2182B",
    midpoint = median(cor_df$rho, na.rm = TRUE),
    name     = "Spearman\ncorrelation"
  ) +
  facet_grid(. ~ species, scales = "free_x", space = "free_x") +
  labs(
    title    = "Bulk RNA-seq vs scRNA-seq pseudobulk: Spearman correlation",
    subtitle = paste(subtitle_parts, collapse = "  |  "),
    x        = "scRNA-seq timepoint",
    y        = "Bulk condition"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 10),
    strip.text       = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    legend.position  = "right"
  )

# Width scales with number of timepoint columns
n_cols <- ncol(cor_all)
pdf(file.path(RESULTS_DIR, "bulk_scrna_correlation_heatmap.pdf"),
    width = max(10, n_cols * 1.1 + 3), height = 5)
print(p)
dev.off()
message("Saved bulk_scrna_correlation_heatmap.pdf")

# =============================================================================
# 7. BEST-MATCH SUMMARY
# =============================================================================

message("\n--- Best-matching scRNA-seq timepoint per bulk condition ---")
for (cond in BULK_ORDER) {
  if (!cond %in% rownames(cor_all)) next
  row  <- cor_all[cond, ]
  best <- names(which.max(row))
  message(sprintf("  %-15s -> %-20s  (rho = %.3f)", cond, best, row[best]))
}
