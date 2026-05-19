# =============================================================================
# project_kawasumi_mds.R
#
# Loads Salmon quantification output from Kawasumi-Kita et al. 2024
# (PRJNA785721), converts to HGNC symbols, averages replicates per condition,
# and projects into the existing no-tadpole centered Spearman MDS space.
#
# Prerequisite: run scripts/kawasumi2024_quantify.sh first to produce
#   data/xenopus_kawasumi2024/quant/<SRR>/quant.sf for all 24 SRRs.
#
# Output:
#   figures/pca_noxentad/kawasumi_projected_PC3vPC1.png
#   results/kawasumi_mds_coords.csv
#
# Run from comparison/ directory: source("R/project_kawasumi_mds.R")
# =============================================================================

suppressPackageStartupMessages({
  library(tximport)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})
source("R/utils_harmonize.R")

CACHE     <- "results/pseudobulk_cache.rds"
QUANT_DIR <- "data/xenopus_kawasumi2024/quant"
TX2GENE   <- "data/xenopus_kawasumi2024/ref/tx2gene.tsv"
FIG_DIR   <- "figures/pca_noxentad"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Check prerequisites
# =============================================================================

if (!file.exists(CACHE))
  stop("Cache missing. Run bulk_vs_scrna_correlation.R first.")
if (!file.exists(TX2GENE))
  stop("tx2gene.tsv missing. Run scripts/kawasumi2024_quantify.sh first.")

srr_dirs <- list.dirs(QUANT_DIR, full.names = TRUE, recursive = FALSE)
quant_files <- file.path(srr_dirs, "quant.sf")
present <- file.exists(quant_files)
if (sum(present) == 0)
  stop("No quant.sf files found in ", QUANT_DIR, ". Run kawasumi2024_quantify.sh first.")
if (sum(!present) > 0)
  message("Warning: ", sum(!present), " quant.sf files missing — proceeding with ",
          sum(present), " of ", length(quant_files))

srr_ids    <- basename(srr_dirs)[present]
quant_files <- quant_files[present]
names(quant_files) <- srr_ids
message("Found ", length(quant_files), " quant.sf files")

# =============================================================================
# 2. Sample manifest — SRR -> condition mapping
# =============================================================================

manifest <- data.frame(
  srr = c(
    "SRR27895473","SRR27895474","SRR27895475","SRR27895476",
    "SRR27908757","SRR27908758","SRR27908759","SRR27908760",
    "SRR27911609","SRR27911610","SRR27911611","SRR27911612",
    "SRR17108583","SRR17108584","SRR17108585","SRR17108587","SRR17108588","SRR17108589",
    "SRR17108577","SRR17108578","SRR17108579","SRR17108580","SRR17108581","SRR17108582"
  ),
  condition = c(
    rep("Kaw_Control_7dpa",  4),
    rep("Kaw_Control_14dpa", 4),
    rep("Kaw_Control_21dpa", 4),
    rep("Kaw_Regen_5dpa",    6),
    rep("Kaw_Regen_7dpa",    6)
  ),
  group = c(
    rep("Xenopus froglet control", 12),
    rep("Xenopus tadpole regen",    12)
  ),
  stringsAsFactors = FALSE
)

# Keep only SRRs with quant files
manifest <- manifest[manifest$srr %in% srr_ids, ]

# =============================================================================
# 3. tximport — transcript -> gene level, then HGNC conversion
# =============================================================================

message("Reading tx2gene table...")
t2g <- read.table(TX2GENE, header = FALSE, sep = "\t",
                   col.names = c("tx_id", "gene_id", "gene_symbol"),
                   stringsAsFactors = FALSE)

# Use gene_symbol column (col1a1.L / col1a1.S format)
t2g_slim <- t2g[, c("tx_id", "gene_symbol")]

message("Running tximport...")
txi <- tximport(
  files      = quant_files[manifest$srr],
  type       = "salmon",
  tx2gene    = t2g_slim,
  ignoreTxVersion = TRUE,
  countsFromAbundance = "lengthScaledTPM"
)

# txi$counts is genes x samples; use TPM-scaled counts for pseudobulk
expr_raw <- txi$counts   # gene rows have X. laevis symbol names (col1a1.L etc.)
colnames(expr_raw) <- manifest$srr

message("Raw expression matrix: ", nrow(expr_raw), " genes x ", ncol(expr_raw), " samples")

# =============================================================================
# 4. Convert X. laevis symbols -> HGNC, aggregate homeologs
# =============================================================================

message("Converting to HGNC symbols...")
expr_hgnc <- aggregate_to_hgnc(expr_raw, species = "xenopus")
message("After HGNC aggregation: ", nrow(expr_hgnc), " genes")

# =============================================================================
# 5. Average replicates per condition -> pseudobulk
# =============================================================================

conditions <- unique(manifest$condition)
pseudo_list <- lapply(conditions, function(cond) {
  srrs <- manifest$srr[manifest$condition == cond]
  srrs <- intersect(srrs, colnames(expr_hgnc))
  mat  <- expr_hgnc[, srrs, drop = FALSE]
  rowMeans(mat)
})
names(pseudo_list) <- conditions

kaw_pseudo <- do.call(cbind, pseudo_list)
message("Pseudobulk matrix: ", nrow(kaw_pseudo), " genes x ", ncol(kaw_pseudo), " conditions")

# =============================================================================
# 6. Load reference MDS — no-tadpole centered space
# =============================================================================

cache   <- readRDS(CACHE)
common  <- cache$common_genes
pseudos <- cache$pseudos

axo_mat <- pseudos[["Axolotl (collaborator)"]][common, ]
mus_mat <- pseudos[["Mouse digit tip"]][common, ]
xen_all <- pseudos[["Xenopus (Lin 2021)"]]
xen_keep <- grepl("^LinBL_", colnames(xen_all)) &
            !colnames(xen_all) %in% c("LinBL_7-14dpa", "LinBL_14-52dpa")
xen_mat <- xen_all[common, xen_keep, drop = FALSE]

ref_mat <- cbind(axo_mat, xen_mat, mus_mat)

ref_species <- c(
  setNames(rep("Axolotl", ncol(axo_mat)), colnames(axo_mat)),
  setNames(rep("Xenopus", ncol(xen_mat)), colnames(xen_mat)),
  setNames(rep("Mouse",   ncol(mus_mat)), colnames(mus_mat))
)

# Within-species centering
ref_centered <- ref_mat
for (sp in unique(ref_species)) {
  cols <- names(ref_species)[ref_species == sp]
  ref_centered[, cols] <- ref_mat[, cols] - rowMeans(ref_mat[, cols, drop = FALSE])
}

# Spearman MDS
rho_ref <- cor(ref_centered, method = "spearman")
mds_fit <- cmdscale(as.dist(1 - rho_ref), k = 6, eig = TRUE)
eig_pos <- pmax(mds_fit$eig, 0)
pct_mds <- 100 * eig_pos[1:6] / sum(eig_pos)

message(sprintf("MDS variance: PC1=%.1f%% PC2=%.1f%% PC3=%.1f%%",
                pct_mds[1], pct_mds[2], pct_mds[3]))

# =============================================================================
# 7. Spearman-weighted centroid projection
# =============================================================================

# Reference Spearman correlation matrix (genes x ref samples, centered)
ref_c_common <- ref_centered  # already on common genes

project_sample <- function(expr_vec) {
  # expr_vec: named numeric vector over HGNC genes
  shared <- intersect(names(expr_vec), rownames(ref_c_common))
  if (length(shared) < 500) {
    warning("Only ", length(shared), " shared genes for projection")
  }
  rho <- cor(expr_vec[shared],
             ref_c_common[shared, ],
             method = "spearman")[1, ]
  # Weighted centroid in 6D MDS space
  w   <- pmax(rho, 0)
  if (sum(w) == 0) w <- rep(1, length(w))
  colSums(mds_fit$points * w) / sum(w)
}

message("Projecting Kawasumi samples...")
kaw_shared <- intersect(rownames(kaw_pseudo), common)
message("  Shared genes with reference: ", length(kaw_shared))

kaw_6d <- t(sapply(conditions, function(cond) {
  project_sample(kaw_pseudo[kaw_shared, cond])
}))
colnames(kaw_6d) <- paste0("PC", 1:6)
message("Projection complete.")

# =============================================================================
# 8. Orient axes (PC1 flipped so high = more regenerative, PC3 = proliferation)
# =============================================================================

ref_coords <- as.data.frame(mds_fit$points)
colnames(ref_coords) <- paste0("PC", 1:6)

# Flip PC1 so regenerating samples score positive
ref_coords$PC1 <- -ref_coords$PC1
kaw_6d[, 1]    <- -kaw_6d[, 1]

# =============================================================================
# 9. Build plot data
# =============================================================================

# Reference sample metadata
ref_labels <- c(
  setNames(
    gsub("AxoCollab_", "", gsub("_", " ", colnames(axo_mat))),
    colnames(axo_mat)
  ),
  setNames(
    gsub("LinBL_", "BL ", colnames(xen_mat)),
    colnames(xen_mat)
  ),
  setNames(
    gsub("Mouse_", "Mouse ", colnames(mus_mat)),
    colnames(mus_mat)
  )
)

ref_groups <- dplyr::case_when(
  grepl("NonRegen", rownames(ref_coords)) ~ "Axolotl (non-regen)",
  grepl("AxoCollab", rownames(ref_coords)) ~ "Axolotl (regenerating)",
  grepl("^LinBL_",  rownames(ref_coords)) ~ "Xenopus froglet (Lin 2021)",
  grepl("^Mouse_",  rownames(ref_coords)) ~ "Mouse digit tip"
)

ref_df <- ref_coords[, c("PC1", "PC3")]
ref_df$label  <- ref_labels[rownames(ref_df)]
ref_df$group  <- ref_groups
ref_df$source <- "Reference"

# Kawasumi projected samples
kaw_labels <- c(
  Kaw_Control_7dpa  = "Control 7dpa",
  Kaw_Control_14dpa = "Control 14dpa",
  Kaw_Control_21dpa = "Control 21dpa",
  Kaw_Regen_5dpa    = "Regen 5dpa",
  Kaw_Regen_7dpa    = "Regen 7dpa"
)

kaw_group <- c(
  Kaw_Control_7dpa  = "Xenopus froglet (Kaw 2024 control)",
  Kaw_Control_14dpa = "Xenopus froglet (Kaw 2024 control)",
  Kaw_Control_21dpa = "Xenopus froglet (Kaw 2024 control)",
  Kaw_Regen_5dpa    = "Xenopus tadpole (Kaw 2024 regen)",
  Kaw_Regen_7dpa    = "Xenopus tadpole (Kaw 2024 regen)"
)

kaw_df <- as.data.frame(kaw_6d[, c("PC1", "PC3")])
kaw_df$label  <- kaw_labels[rownames(kaw_df)]
kaw_df$group  <- kaw_group[rownames(kaw_df)]
kaw_df$source <- "Kawasumi 2024"

plot_df <- rbind(ref_df, kaw_df)

# =============================================================================
# 10. Plot
# =============================================================================

group_colors <- c(
  "Axolotl (regenerating)"             = "#1565C0",
  "Axolotl (non-regen)"                = "#90CAF9",
  "Xenopus froglet (Lin 2021)"         = "#F4511E",
  "Mouse digit tip"                    = "#43A047",
  "Xenopus froglet (Kaw 2024 control)" = "#E65100",
  "Xenopus tadpole (Kaw 2024 regen)"   = "#00695C"
)

group_shapes <- c(
  "Axolotl (regenerating)"             = 21,
  "Axolotl (non-regen)"                = 24,
  "Xenopus froglet (Lin 2021)"         = 21,
  "Mouse digit tip"                    = 21,
  "Xenopus froglet (Kaw 2024 control)" = 23,
  "Xenopus tadpole (Kaw 2024 regen)"   = 23
)

p <- ggplot(plot_df, aes(PC3, PC1, fill = group, shape = group, label = label)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
  geom_point(data = subset(plot_df, source == "Reference"),
             size = 5, color = "white", stroke = 0.5, alpha = 0.95) +
  geom_point(data = subset(plot_df, source == "Kawasumi 2024"),
             size = 5, color = "black", stroke = 0.8, alpha = 0.90) +
  geom_text_repel(size = 3.2, color = "grey20", bg.color = "white", bg.r = 0.08,
                  box.padding = 0.5, point.padding = 0.3,
                  max.overlaps = 30,
                  segment.color = "grey60", segment.size = 0.3) +
  scale_fill_manual(values  = group_colors, name = NULL) +
  scale_shape_manual(values = group_shapes, name = NULL) +
  labs(
    title    = "Cross-species regeneration MDS with Kawasumi 2024 projection",
    subtitle = sprintf(
      "Spearman MDS (no tadpole reference) | PC1=%.1f%% PC3=%.1f%% | diamonds = Kawasumi 2024",
      pct_mds[1], pct_mds[3]
    ),
    x = sprintf("PC3 - proliferation/fibrosis axis (%.1f%%)", pct_mds[3]),
    y = sprintf("PC1 - regenerative gradient (%.1f%%)", pct_mds[1])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 14, color = "grey15"),
    plot.subtitle     = element_text(size = 9, color = "grey45"),
    axis.title        = element_text(size = 11, color = "grey25"),
    axis.text         = element_text(size = 9, color = "grey40"),
    legend.text       = element_text(size = 10, color = "grey20"),
    legend.position   = "right",
    legend.key        = element_rect(fill = NA, color = NA),
    plot.margin       = margin(16, 16, 16, 16)
  ) +
  guides(fill  = guide_legend(override.aes = list(size = 4.5)),
         shape = guide_legend(override.aes = list(size = 4.5)))

out_png <- file.path(FIG_DIR, "kawasumi_projected_PC3vPC1.png")
ggsave(out_png, p, width = 11, height = 8, dpi = 200, bg = "white")
message("Saved: ", out_png)

# =============================================================================
# 11. Save coordinates
# =============================================================================

all_coords <- rbind(
  cbind(ref_df, dataset = "Reference"),
  cbind(kaw_df, dataset = "Kawasumi2024")
)
out_csv <- "results/kawasumi_mds_coords.csv"
write.csv(all_coords, out_csv, row.names = TRUE)
message("Saved: ", out_csv)

cat("\nKawasumi 2024 projected positions:\n")
print(round(kaw_df[, c("PC3","PC1")], 3))
