# =============================================================================
# 02 — Trinity ID annotation, feature renaming, Harmony integration
# Input:  results/axo_{intact,wh,eb,mb}.rds
# Output: results/axo_integrated.rds
#         helper/trinity_to_symbol_map.csv
#
# Trinity IDs embed UniProt annotations: c#####_g#_i#^sp|ACCESSION|SYMBOL_SPECIES^
# Seurat sanitizes | and _ to -, so regex patterns use - throughout.
# For GOI with non-obvious UniProt entry names, UNIPROT_ALIASES provides fallbacks.
# =============================================================================

library(Seurat)
library(harmony)
library(dplyr)
library(stringr)
library(Matrix)

results_dir <- "results"
helper_dir  <- "helper"
dir.create(helper_dir, showWarnings = FALSE)

TIMEPOINT_ORDER <- c("Intact", "WH_N4", "WH_N5", "WH_N6", "EB", "MB")

GOI <- c(
  "GLI1","PTCH1","PTCH2","SHH","IHH","DHH",
  "GDF5","MSX1","SALL1","GREM1",
  "FREM1","CYTL1","CLDN5","IL11","FGFBP1",
  "C3","LDLR","ALOX15B","MSMO1","CYP51A1","FDPS"
)

# UniProt entry-name aliases for GOI that don't match HGNC symbol directly
UNIPROT_ALIASES <- list(
  "PTCH1"   = c("PTCH1","PTC1"),
  "PTCH2"   = c("PTCH2","PTC2"),
  "SALL1"   = c("SALL1","SAL1"),
  "CLDN5"   = c("CLDN5","CLD5"),
  "FGFBP1"  = c("FGFBP1","FGBP1"),
  "C3"      = c("C3","CO3"),
  "ALOX15B" = c("ALOX15B","LX15B"),
  "CYP51A1" = c("CYP51A1","CP51A"),
  "GLI1"    = c("GLI1"),
  "SHH"     = c("SHH"),
  "IHH"     = c("IHH"),
  "DHH"     = c("DHH"),
  "GDF5"    = c("GDF5"),
  "MSX1"    = c("MSX1","HME1"),
  "GREM1"   = c("GREM1"),
  "FREM1"   = c("FREM1"),
  "CYTL1"   = c("CYTL1"),
  "IL11"    = c("IL11"),
  "LDLR"    = c("LDLR","LDLR1"),
  "MSMO1"   = c("MSMO1"),
  "FDPS"    = c("FDPS")
)

# =============================================================================
# 1. LOAD OBJECTS
# =============================================================================

intact <- readRDS(file.path(results_dir, "axo_intact.rds"))
wh     <- readRDS(file.path(results_dir, "axo_wh.rds"))
eb     <- readRDS(file.path(results_dir, "axo_eb.rds"))
mb     <- readRDS(file.path(results_dir, "axo_mb.rds"))

# =============================================================================
# 2. BUILD TRINITY -> GENE SYMBOL MAP
# =============================================================================

message("\nBuilding Trinity ID -> gene symbol map...")
all_genes <- Reduce(union, list(
  rownames(intact), rownames(wh), rownames(eb), rownames(mb)
))
message("Total unique genes across all studies: ", length(all_genes))

# Extract embedded UniProt gene symbol from sanitized Trinity ID
# Sanitized format: ...^sp-ACCESSION-GENESYM-SPECIES^...
extract_embedded_symbol <- function(ids) {
  str_extract(ids, "(?<=-)[A-Z][A-Z0-9]{1,}(?=-[A-Z]{2,5}\\^)")
}

embedded_syms <- extract_embedded_symbol(all_genes)
embedded_map <- data.frame(
  trinity_id  = all_genes,
  gene_symbol = embedded_syms,
  source      = "embedded",
  stringsAsFactors = FALSE
) |> filter(!is.na(gene_symbol))

message("  Embedded annotations: ", nrow(embedded_map))

# GOI-specific direct search using UNIPROT_ALIASES
goi_direct_map <- bind_rows(lapply(names(UNIPROT_ALIASES), function(hgnc) {
  aliases <- UNIPROT_ALIASES[[hgnc]]
  pattern <- paste0("-(?:", paste(aliases, collapse = "|"), ")-[A-Z]{2,5}\\^")
  hits    <- grep(pattern, all_genes, value = TRUE, ignore.case = TRUE)
  if (length(hits) > 0)
    data.frame(trinity_id  = hits,
               gene_symbol = hgnc,
               source      = "goi_direct",
               stringsAsFactors = FALSE)
}))
message("  GOI direct matches: ", nrow(goi_direct_map),
        " Trinity IDs covering ",
        length(unique(goi_direct_map$gene_symbol)), " GOI")

# Merge: goi_direct overrides embedded for GOI, embedded for everything else
gene_map <- bind_rows(goi_direct_map, embedded_map) |>
  filter(trinity_id %in% all_genes) |>
  mutate(source_rank = if_else(source == "goi_direct", 1L, 2L)) |>
  arrange(trinity_id, source_rank) |>
  group_by(trinity_id) |>
  slice(1) |>
  ungroup()

message("Total mappings: ", nrow(gene_map))

goi_coverage <- gene_map |> filter(gene_symbol %in% GOI)
message("GOI found: ", length(unique(goi_coverage$gene_symbol)), "/", length(GOI))

write.csv(gene_map, file.path(helper_dir, "trinity_to_symbol_map.csv"),
          row.names = FALSE)
message("Saved helper/trinity_to_symbol_map.csv")

# =============================================================================
# 3. RENAME FEATURES & NORMALISE (per dataset)
# =============================================================================

rename_and_normalize <- function(obj, gene_map) {
  current_genes <- rownames(obj)
  rename_vec    <- setNames(gene_map$gene_symbol, gene_map$trinity_id)
  rename_vec    <- rename_vec[names(rename_vec) %in% current_genes]

  # Skip renaming if new name already exists under a different Trinity ID
  already_present <- rename_vec %in% current_genes & names(rename_vec) != rename_vec
  rename_vec <- rename_vec[!already_present]

  new_names <- current_genes
  idx <- match(names(rename_vec), current_genes)
  new_names[idx[!is.na(idx)]] <- rename_vec[!is.na(idx)]

  # Deduplicate: when multiple Trinity IDs map to the same symbol,
  # keep the one with highest total counts
  dup_syms <- unique(new_names[duplicated(new_names)])
  if (length(dup_syms) > 0) {
    counts_tmp <- GetAssayData(obj, layer = "counts")
    for (sym in dup_syms) {
      dup_idx    <- which(new_names == sym)
      row_totals <- Matrix::rowSums(counts_tmp[dup_idx, ])
      keep       <- dup_idx[which.max(row_totals)]
      revert     <- dup_idx[dup_idx != keep]
      new_names[revert] <- current_genes[revert]
    }
  }

  counts_mat <- GetAssayData(obj, layer = "counts")
  rownames(counts_mat) <- new_names
  obj[["RNA"]] <- CreateAssay5Object(counts = counts_mat)

  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj
}

message("\nRenaming features and normalizing...")
intact <- rename_and_normalize(intact, gene_map)
wh     <- rename_and_normalize(wh,     gene_map)
eb     <- rename_and_normalize(eb,     gene_map)
mb     <- rename_and_normalize(mb,     gene_map)

# =============================================================================
# 4. MERGE & HARMONY INTEGRATION
# =============================================================================

message("\nMerging all objects...")
axo_merged <- merge(
  intact,
  y            = list(wh, eb, mb),
  add.cell.ids = c("Intact", "WH", "EB", "MB"),
  merge.data   = TRUE
)
axo_merged$timepoint <- factor(axo_merged$timepoint, levels = TIMEPOINT_ORDER)

message("Merged: ", ncol(axo_merged), " cells, ", nrow(axo_merged), " genes")
print(table(axo_merged$timepoint))

axo_merged <- FindVariableFeatures(axo_merged, nfeatures = 2000, verbose = FALSE)
axo_merged <- ScaleData(axo_merged,
                         features = VariableFeatures(axo_merged), verbose = FALSE)
axo_merged <- RunPCA(axo_merged, npcs = 50, verbose = FALSE)

message("Running Harmony...")
axo_merged <- RunHarmony(
  axo_merged,
  group.by.vars  = "study_id",
  reduction      = "pca",
  reduction.save = "harmony",
  verbose        = FALSE
)

message("Computing integrated UMAP...")
axo_merged <- RunUMAP(axo_merged, reduction = "harmony", dims = 1:30,
                       reduction.name = "umap_harmony", verbose = FALSE)

axo_merged <- FindNeighbors(axo_merged, reduction = "harmony", dims = 1:30,
                              verbose = FALSE)
axo_merged <- FindClusters(axo_merged, resolution = 0.5, verbose = FALSE)

# =============================================================================
# 5. SAVE
# =============================================================================

DefaultDimReduc(axo_merged) <- "umap_harmony"
saveRDS(axo_merged, file.path(results_dir, "axo_integrated.rds"))
message("Saved results/axo_integrated.rds")

message("\nGOI check in integrated object:")
goi_present <- intersect(GOI, rownames(axo_merged))
message("  Found (", length(goi_present), "/", length(GOI), "): ",
        paste(sort(goi_present), collapse = ", "))
message("  Missing: ", paste(setdiff(GOI, goi_present), collapse = ", "))
