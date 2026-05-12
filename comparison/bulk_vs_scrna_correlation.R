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
XEN_RDS   <- "~/Desktop/scrnaseq-pipelines/xenopus/results/xen_merged.rds"
MOUSE_RDS        <- "~/Desktop/scrnaseq-pipelines/mouse-digit/results/mouse_merged.rds"
MOUSE_NONREGEN_RDS <- "~/Desktop/scrnaseq-pipelines/mouse-digit-nonregen/results/mouse_nonregen_merged.rds"
AXO_COLLAB_DIR     <- "data/axo_collab"   # ShinyCell export from collaborator dataset

RESULTS_DIR <- "results"
dir.create(RESULTS_DIR, showWarnings = FALSE)

BULK_ORDER  <- c("BMS_acute",   "BMS_24hrs",   "Pretreatment",   "Reamputation",
                 "EtOH_acute",  "EtOH_24hrs",  "EtOH_pretreat",  "EtOH_reamp")

# =============================================================================
# 1. BULK DATA
# =============================================================================

message("Loading bulk RNA-seq data...")

# col_pattern: grep pattern to select sample columns (e.g. "BMS" or "EtOH").
# NULL averages all columns (original behaviour, but should not be used now
# that each file contains both BMS and EtOH replicates).
.load_bulk_condition <- function(path, condition, col_pattern = NULL) {
  df    <- read.csv(path, row.names = 1, check.names = FALSE)
  named <- df[!grepl("^LOC", rownames(df)), , drop = FALSE]
  if (!is.null(col_pattern))
    named <- named[, grepl(col_pattern, colnames(named), ignore.case = TRUE), drop = FALSE]
  setNames(data.frame(rowMeans(named)), condition)
}

bulk_list <- list(
  # BMS-treated replicates only
  .load_bulk_condition(file.path(BULK_DIR, "BMS_normalized_counts.csv"),          "BMS_acute",    "BMS"),
  .load_bulk_condition(file.path(BULK_DIR, "24 hrs_normalized_counts.csv"),        "BMS_24hrs",    "BMS"),
  .load_bulk_condition(file.path(BULK_DIR, "Pretreatment_normalized_counts.csv"),  "Pretreatment", "BMS"),
  .load_bulk_condition(file.path(BULK_DIR, "Reamputation_normalized_counts.csv"),  "Reamputation", "BMS"),
  # EtOH vehicle controls (paired with BMS conditions above)
  .load_bulk_condition(file.path(BULK_DIR, "BMS_normalized_counts.csv"),          "EtOH_acute",   "EtOH"),
  .load_bulk_condition(file.path(BULK_DIR, "24 hrs_normalized_counts.csv"),        "EtOH_24hrs",   "EtOH"),
  .load_bulk_condition(file.path(BULK_DIR, "Pretreatment_normalized_counts.csv"),  "EtOH_pretreat","EtOH"),
  .load_bulk_condition(file.path(BULK_DIR, "Reamputation_normalized_counts.csv"),  "EtOH_reamp",   "EtOH")
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

# Strip heavy layers and join before pseudobulk extraction (halves peak RAM).
.slim_for_pseudobulk <- function(obj) {
  for (l in Layers(obj)[startsWith(Layers(obj), "counts.")])
    LayerData(obj, layer = l) <- NULL
  for (l in Layers(obj)[startsWith(Layers(obj), "scale.data")])
    LayerData(obj, layer = l) <- NULL
  gc()
  JoinLayers(obj)
}

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

# --- Xenopus laevis Lin 2021 (pseudotetraploid, .L/.S homeologs) ------------
#
# The object contains two experimental groups (set in 01_xenopus_qc_clustering.R):
#   condition == "BL"   : froglet post-amputation, non-regenerative (dpa labels)
#   condition == "LBst" : tadpole limb bud stage, regenerative (NF stage labels)
# Each group gets its own column prefix so they appear as separate heatmap facets:
#   LinBL_   -> "Xenopus froglet (Lin 2021)"
#   LinLBst_ -> "Xenopus tadpole (Lin 2021)"

pseudobulk_xenopus_lin <- function(rds_path) {
  message("Computing pseudobulk: Xenopus (Lin 2021)...")
  xen  <- .slim_for_pseudobulk(readRDS(rds_path))
  expr <- GetAssayData(xen, layer = "data")

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

  # --- Froglet BL (non-regenerative) ---
  bl_idx    <- which(xen$condition == "BL")
  bl_order  <- c("0dpa","3dpa","7-14dpa","14dpa","14-52dpa")
  bl_factor <- factor(xen$timepoint[bl_idx], levels = bl_order)
  pb_bl <- .mean_per_tp(expr_summed[, bl_idx], bl_factor, bl_order, prefix = "LinBL")

  # --- Tadpole LBst (regenerative) ---
  lb_idx    <- which(xen$condition == "LBst")
  lb_order  <- c("NF50","NF51","NF52")
  lb_factor <- factor(xen$timepoint[lb_idx], levels = lb_order)
  pb_lb <- .mean_per_tp(expr_summed[, lb_idx], lb_factor, lb_order, prefix = "LinLBst")

  message("  BL cells: ", length(bl_idx), " | LBst cells: ", length(lb_idx))
  as.matrix(cbind(pb_bl, pb_lb))
}

# --- Mouse digit tip regenerative (Johnson, Masias & Lehoczky 2020, GSE143888) ---

pseudobulk_mouse_digit <- function(rds_path) {
  message("Computing pseudobulk: Mouse digit tip (Johnson/Lehoczky 2020, GSE143888)...")
  mus       <- .slim_for_pseudobulk(readRDS(rds_path))
  expr      <- GetAssayData(mus, layer = "data")
  tp_order  <- c("0dpa","11dpa","12dpa","14dpa","17dpa")
  tp_factor <- factor(mus$timepoint, levels = tp_order)

  pb <- .mean_per_tp(expr, tp_factor, tp_order, prefix = "Mouse")
  message("  Mouse: ", ncol(expr), " cells x ", nrow(expr), " genes")

  aggregate_to_hgnc(as.matrix(pb), species = "mouse")
}

# --- Collaborator axolotl dataset (private/unpublished) ---------------------
#
# Data stored as ShinyCell export in data/axo_collab/:
#   merged_objmeta.rds  — cell metadata (43,083 cells)
#   merged_objgene.rds  — named integer vector: HGNC symbol -> h5 row index
#   merged_objassay_RNA.h5 — dense float matrix: genes x cells (grp/data)
#
# Conditions:
#   Regenerating   : dpa3 (wound healing), dpa14 (early bud), dpa23 (medium bud)
#   Non-regenerating: labeled "limb" in timepoint column (proximal amputation controls)
#
# Genes are already HGNC symbols — no Trinity ID mapping required.

pseudobulk_axo_collab <- function(data_dir, chunk_size = 1000) {
  message("Computing pseudobulk: Collaborator axolotl...")
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r package required: install.packages('hdf5r')")

  meta     <- readRDS(file.path(data_dir, "merged_objmeta.rds"))
  gene_idx <- readRDS(file.path(data_dir, "merged_objgene.rds"))[["RNA"]]

  # Gene names in h5 row order (values are row positions 1..n)
  gene_names <- names(sort(gene_idx))
  n_genes    <- length(gene_names)

  # Derive condition from sample_id
  meta$condition <- ifelse(grepl("^non-regenerating", meta$sample_id),
                           "NonRegen", "Regen")
  meta$group <- paste0(meta$condition, "_", meta$timepoint)

  groups <- c("Regen_dpa3", "Regen_dpa14", "Regen_dpa23", "NonRegen_limb")
  group_cell_idx <- lapply(groups, function(g) which(meta$group == g))
  names(group_cell_idx) <- groups

  message("  Cells per group: ",
          paste(names(group_cell_idx),
                sapply(group_cell_idx, length), sep = "=", collapse = ", "))

  # Accumulate pseudobulk means reading h5 in gene chunks
  pb <- matrix(0, nrow = n_genes, ncol = length(groups),
               dimnames = list(gene_names, groups))

  h5   <- hdf5r::H5File$new(file.path(data_dir, "merged_objassay_RNA.h5"), mode = "r")
  dset <- h5[["grp/data"]]
  n_chunks <- ceiling(n_genes / chunk_size)

  for (i in seq_len(n_chunks)) {
    g_start <- (i - 1L) * chunk_size + 1L
    g_end   <- min(i * chunk_size, n_genes)
    chunk   <- dset[g_start:g_end, ]   # (chunk_size x n_cells) dense matrix
    for (g in groups) {
      idx <- group_cell_idx[[g]]
      if (length(idx) > 0L)
        pb[g_start:g_end, g] <- rowMeans(chunk[, idx, drop = FALSE])
    }
    if (i %% 10 == 0)
      message("  chunk ", i, "/", n_chunks)
  }
  h5$close_all()

  message("  Done: ", n_genes, " HGNC genes x ", length(groups), " groups")
  colnames(pb) <- paste0("AxoCollab_", colnames(pb))
  pb
}

# --- Mouse digit: regenerative vs non-regenerative (Storer 2020, GSE135985) -
#
# Same lab as GSE143888. Object contains both conditions (condition column =
# "regenerative" | "non-regenerative"). We build two separate pseudobulk
# columns prefixed by "MouseRegen_" and "MouseNonRegen_" for direct comparison.

pseudobulk_mouse_nonregen <- function(rds_path) {
  if (!file.exists(path.expand(rds_path))) {
    message("  SKIPPING mouse non-regen: RDS not found (run mouse-digit-nonregen pipeline first)")
    return(NULL)
  }
  message("Computing pseudobulk: Mouse digit non-regen (Storer 2020, GSE135985)...")
  mus  <- .slim_for_pseudobulk(readRDS(rds_path))
  expr <- GetAssayData(mus, layer = "data")
  message("  Mouse non-regen: ", ncol(expr), " cells x ", nrow(expr), " genes")

  regen_idx    <- which(mus$condition == "regenerative")
  nonregen_idx <- which(mus$condition == "non-regenerative")
  tp_regen     <- c("uninjured","7dpa","10dpa","14dpa","28dpa","56dpa")
  tp_nonregen  <- c("10dpa","14dpa")

  pb_regen    <- .mean_per_tp(expr[, regen_idx],
                               factor(mus$timepoint[regen_idx],    levels = tp_regen),
                               tp_regen,    prefix = "MouseRegen")
  pb_nonregen <- .mean_per_tp(expr[, nonregen_idx],
                               factor(mus$timepoint[nonregen_idx], levels = tp_nonregen),
                               tp_nonregen, prefix = "MouseNonRegen")
  pb <- cbind(pb_regen, pb_nonregen)
  aggregate_to_hgnc(as.matrix(pb), species = "mouse")
}

# =============================================================================
# 4. BUILD PSEUDOBULK LIST
#    Each entry: HGNC-gene-symbol matrix x timepoints
# =============================================================================

mouse_nonregen_pb <- pseudobulk_mouse_nonregen(MOUSE_NONREGEN_RDS)

pseudos <- Filter(Negate(is.null), list(
  "Axolotl (collaborator)"      = pseudobulk_axo_collab(AXO_COLLAB_DIR),
  "Xenopus (Lin 2021)"          = pseudobulk_xenopus_lin(XEN_RDS),
  "Mouse digit tip"             = pseudobulk_mouse_digit(MOUSE_RDS),
  "Mouse regen/non-regen"       = mouse_nonregen_pb
))

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

message("\nCorrelation matrix (per-pairwise shared genes):")
print(round(cor_all, 3))

write.csv(round(cor_all, 4),
          file.path(RESULTS_DIR, "bulk_scrna_correlation_matrix.csv"))
message("Saved bulk_scrna_correlation_matrix.csv")

# =============================================================================
# 5b. COMMON GENE SET CORRELATION
#
# Re-run Spearman using only genes shared across ALL datasets simultaneously
# (bulk ∩ axolotl ∩ xenopus ∩ mouse). This removes the gene-count advantage
# that Xenopus/Mouse have over Axolotl (5,617 vs ~10,000 genes), making rho
# values directly comparable across species.
# =============================================================================

common_genes <- Reduce(intersect, c(
  list(rownames(bulk_mat)),
  lapply(pseudos, rownames)
))
message("\nCommon gene set (bulk + all scRNA-seq datasets): ", length(common_genes), " genes")

bulk_common <- bulk_mat[common_genes, , drop = FALSE]
cor_list_common <- lapply(names(pseudos), function(nm) {
  s <- pseudos[[nm]][common_genes, , drop = FALSE]
  cor_out <- matrix(NA_real_, nrow = ncol(bulk_common), ncol = ncol(s),
                    dimnames = list(colnames(bulk_common), colnames(s)))
  for (i in colnames(bulk_common))
    for (j in colnames(s))
      cor_out[i, j] <- cor(bulk_common[, i], s[, j],
                           method = "spearman", use = "complete.obs")
  cor_out
})
names(cor_list_common) <- names(pseudos)
cor_common <- do.call(cbind, cor_list_common)

message("\nCorrelation matrix (common gene set, n = ", length(common_genes), " genes):")
print(round(cor_common, 3))

write.csv(round(cor_common, 4),
          file.path(RESULTS_DIR, "bulk_scrna_correlation_matrix_commongenes.csv"))
message("Saved bulk_scrna_correlation_matrix_commongenes.csv")

# =============================================================================
# 6. HEATMAPS
# =============================================================================

# Derive friendly species labels from column name prefixes
.prefix_to_species <- function(col_names) {
  prefixes <- sub("_.*", "", col_names)
  species_labels <- c(
    AxoCollab      = "Axolotl (collaborator)",
    LinBL          = "Xenopus froglet (Lin 2021)",
    LinLBst        = "Xenopus tadpole (Lin 2021)",
    Mouse          = "Mouse digit tip",
    MouseRegen     = "Mouse regen (GSE135985)",
    MouseNonRegen  = "Mouse non-regen (GSE135985)"
  )
  ifelse(prefixes %in% names(species_labels),
         species_labels[prefixes], prefixes)
}

.make_heatmap <- function(cor_mat, title, subtitle) {
  df <- as.data.frame(cor_mat) |>
    tibble::rownames_to_column("bulk_condition") |>
    pivot_longer(-bulk_condition,
                 names_to  = "scrna_timepoint",
                 values_to = "rho")

  df$bulk_condition  <- factor(df$bulk_condition,  levels = BULK_ORDER)
  df$scrna_timepoint <- factor(df$scrna_timepoint, levels = colnames(cor_mat))
  df$species         <- .prefix_to_species(as.character(df$scrna_timepoint))
  df$species         <- factor(df$species,
                               levels = unique(.prefix_to_species(colnames(cor_mat))))
  # Row group: BMS vs EtOH, for visual separation via row facet
  df$treatment_group <- ifelse(grepl("^EtOH", as.character(df$bulk_condition)),
                               "EtOH (vehicle)", "BMS-345541")
  df$treatment_group <- factor(df$treatment_group,
                               levels = c("BMS-345541", "EtOH (vehicle)"))

  ggplot(df, aes(x = scrna_timepoint, y = bulk_condition, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", rho)), size = 3, color = "black") +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "#F7F7F7",
      high     = "#B2182B",
      midpoint = median(df$rho, na.rm = TRUE),
      name     = "Spearman\ncorrelation"
    ) +
    facet_grid(treatment_group ~ species, scales = "free", space = "free") +
    labs(title = title, subtitle = subtitle,
         x = "scRNA-seq timepoint", y = "Bulk condition") +
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
}

# Page 1: per-pairwise shared genes (original approach)
gene_counts    <- sapply(names(pseudos), function(nm)
  length(intersect(rownames(bulk_mat), rownames(pseudos[[nm]]))))
subtitle_full  <- paste0(names(pseudos), ": ", gene_counts, " genes",
                         collapse = "  |  ")

p_full <- .make_heatmap(
  cor_all,
  title    = "Bulk RNA-seq vs scRNA-seq pseudobulk: Spearman rho (per-pairwise shared genes)",
  subtitle = subtitle_full
)

# Page 2: common gene set — rho values are directly comparable across species
p_common <- .make_heatmap(
  cor_common,
  title    = sprintf("Bulk RNA-seq vs scRNA-seq pseudobulk: Spearman rho (common gene set, n = %d)",
                     length(common_genes)),
  subtitle = "Same genes used for all species — rho directly comparable across Axolotl / Xenopus / Mouse"
)

n_cols <- ncol(cor_all)
n_rows <- nrow(cor_all)
pdf(file.path(RESULTS_DIR, "bulk_scrna_correlation_heatmap.pdf"),
    width = max(10, n_cols * 1.1 + 3), height = max(5, n_rows * 0.6 + 2))
print(p_full)
print(p_common)
dev.off()
message("Saved bulk_scrna_correlation_heatmap.pdf (page 1: per-pairwise; page 2: common gene set)")

# =============================================================================
# 7. BEST-MATCH SUMMARY
# =============================================================================

message("\n--- Best-matching scRNA-seq timepoint per bulk condition (common gene set) ---")
for (cond in BULK_ORDER) {
  if (!cond %in% rownames(cor_common)) next
  row  <- cor_common[cond, ]
  best <- names(which.max(row))
  message(sprintf("  %-15s -> %-20s  (rho = %.3f)", cond, best, row[best]))
}

# Save pseudobulk list so DEG and GOI analysis scripts can load it without
# re-running the expensive computation.
saveRDS(list(pseudos = pseudos, bulk_mat = bulk_mat, common_genes = common_genes),
        file.path(RESULTS_DIR, "pseudobulk_cache.rds"))
message("Saved pseudobulk_cache.rds (load in deg_scrna_analysis.R)")
