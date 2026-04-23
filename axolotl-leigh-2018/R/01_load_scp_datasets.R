# =============================================================================
# 01 — Load Broad Single Cell Portal datasets into Seurat objects
# Input:  data/SCP422_intact/, data/SCP489_WH/, data/SCP499_EB/, data/SCP500_MB/
# Output: results/axo_intact.rds, axo_wh.rds, axo_eb.rds, axo_mb.rds
#
# Source: Leigh et al. 2018, Nat Commun, doi:10.1038/s41467-018-07604-0
#         Broad SCP422/489/499/500 — GEO GSE121737 / SRA SRP167700
#
# SCP matrix format: genes x cells, tab-separated, gzipped
# Idents / coordinates files: 2-row SCP header (NAME row + TYPE row) then data
# =============================================================================

library(Seurat)
library(data.table)
library(stringr)
library(Matrix)

# Run from leigh-2018-axolotl-integration/
data_dir    <- "data"
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# HELPER: load one SCP dataset -> Seurat object
# =============================================================================

load_scp_dataset <- function(matrix_path, idents_path, coords_path,
                              study_id, timepoint) {
  message("\nLoading ", study_id, " (", timepoint, ")...")

  # --- Expression matrix (gzipped tab-delimited genes x cells) ---------------
  message("  Reading matrix...")
  mat_dt <- fread(
    cmd         = paste0("gunzip -c ", shQuote(matrix_path)),
    sep         = "\t",
    header      = TRUE,
    check.names = FALSE,
    data.table  = TRUE
  )
  genes     <- mat_dt[[1]]
  cell_cols <- colnames(mat_dt)[-1]
  mat_dt[, 1 := NULL]

  # Convert in 1,000-cell chunks to limit peak memory
  message("  Converting to sparse matrix (chunked)...")
  chunk_size <- 1000
  chunks     <- split(seq_along(cell_cols),
                       ceiling(seq_along(cell_cols) / chunk_size))
  chunk_mats <- lapply(seq_along(chunks), function(ci) {
    idx   <- chunks[[ci]]
    cols  <- cell_cols[idx]
    chunk <- as.matrix(mat_dt[, ..idx])
    rownames(chunk) <- genes
    colnames(chunk) <- cols
    as(chunk, "dgCMatrix")
  })
  mat <- do.call(cbind, chunk_mats)
  rm(mat_dt, chunk_mats); gc()
  message("  Matrix: ", nrow(mat), " genes x ", ncol(mat), " cells")

  # --- Cell type annotations --------------------------------------------------
  # SCP format: header row (NAME ...) then TYPE row then data rows
  idents_raw <- read.table(idents_path, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)
  idents_raw <- idents_raw[idents_raw[[1]] != "TYPE", , drop = FALSE]
  rownames(idents_raw) <- idents_raw[[1]]
  idents_raw <- idents_raw[, -1, drop = FALSE]

  # --- Published UMAP coordinates ---------------------------------------------
  coords_raw <- read.table(coords_path, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)
  coords_raw <- coords_raw[coords_raw[[1]] != "TYPE", , drop = FALSE]
  rownames(coords_raw) <- coords_raw[[1]]
  coords_raw <- coords_raw[, -1, drop = FALSE]
  coords_raw[[1]] <- as.numeric(coords_raw[[1]])
  coords_raw[[2]] <- as.numeric(coords_raw[[2]])

  # --- Align cells across all three sources -----------------------------------
  shared_cells <- Reduce(intersect, list(
    colnames(mat),
    rownames(idents_raw),
    rownames(coords_raw)
  ))
  message("  Shared cells: ", length(shared_cells))

  mat        <- mat[, shared_cells]
  idents_raw <- idents_raw[shared_cells, , drop = FALSE]
  coords_raw <- coords_raw[shared_cells, , drop = FALSE]

  # --- Build Seurat object ----------------------------------------------------
  obj <- CreateSeuratObject(
    counts       = mat,
    project      = study_id,
    min.cells    = 5,
    min.features = 200
  )

  obj$study_id      <- study_id
  obj$timepoint     <- timepoint
  obj$paper_cluster <- idents_raw[[1]]

  # SCP489 encodes three WH replicates as N4_ / N5_ / N6_ barcode prefixes
  if (study_id == "SCP489") {
    obj$timepoint <- paste0("WH_", str_extract(colnames(obj), "^N[0-9]+"))
  }

  # Embed paper UMAP as a named DimReduc so it can be plotted directly
  paper_umap <- as.matrix(coords_raw)
  colnames(paper_umap) <- c("UMAP_1", "UMAP_2")
  obj[["paper_umap"]] <- CreateDimReducObject(
    embeddings = paper_umap,
    key        = "UMAP_",
    assay      = "RNA"
  )

  message("  Cell types: ",
          paste(sort(unique(obj$paper_cluster)), collapse = ", "))
  obj
}

# =============================================================================
# LOAD FOUR TIMEPOINTS
# =============================================================================

intact <- load_scp_dataset(
  matrix_path = file.path(data_dir, "SCP422_intact/intact.Matrix.txt.gz"),
  idents_path = file.path(data_dir, "SCP422_intact/intact.idents.txt"),
  coords_path = file.path(data_dir, "SCP422_intact/intact.coordinates.txt"),
  study_id    = "SCP422",
  timepoint   = "Intact"
)
gc()

wh <- load_scp_dataset(
  matrix_path = file.path(data_dir, "SCP489_WH/WH.matrix.txt.gz"),
  idents_path = file.path(data_dir, "SCP489_WH/WH.idents.txt"),
  coords_path = file.path(data_dir, "SCP489_WH/WH.coordinates.txt"),
  study_id    = "SCP489",
  timepoint   = "WH"
)
gc()

eb <- load_scp_dataset(
  matrix_path = file.path(data_dir, "SCP499_EB/EB.matrix.txt.gz"),
  idents_path = file.path(data_dir, "SCP499_EB/EB.idents.txt"),
  coords_path = file.path(data_dir, "SCP499_EB/EB.coordinates.txt"),
  study_id    = "SCP499",
  timepoint   = "EB"
)
gc()

mb <- load_scp_dataset(
  matrix_path = file.path(data_dir, "SCP500_MB/MB.matrix.txt.gz"),
  idents_path = file.path(data_dir, "SCP500_MB/MB.idents.txt"),
  coords_path = file.path(data_dir, "SCP500_MB/MB.coordinates.txt"),
  study_id    = "SCP500",
  timepoint   = "MB"
)
gc()

saveRDS(intact, file.path(results_dir, "axo_intact.rds"))
saveRDS(wh,     file.path(results_dir, "axo_wh.rds"))
saveRDS(eb,     file.path(results_dir, "axo_eb.rds"))
saveRDS(mb,     file.path(results_dir, "axo_mb.rds"))
message("\nSaved four Seurat objects to results/")
