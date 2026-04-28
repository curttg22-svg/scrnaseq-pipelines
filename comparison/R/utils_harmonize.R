# =============================================================================
# utils_harmonize.R — Gene symbol harmonization to HGNC uppercase
#
# Purpose: convert species-specific gene naming conventions to a common
#   HGNC uppercase namespace so scRNA-seq datasets can be compared against
#   axolotl bulk RNA-seq (which uses current NCBI/HGNC symbols).
#
# Supported species:
#   "human"   — NCBI/HGNC symbols (the bulk format). Passes each name
#               through org.Hs.eg.db's ALIAS table, which covers official
#               symbols, NCBI-listed synonyms, and some SwissProt entry
#               names recorded as gene aliases. Catches partial improvements
#               for Leigh 2018 axolotl Trinity data (Trinity IDs annotated
#               against SwissProt; some abbreviations differ from HGNC).
#               Falls back to toupper() for names not in the alias table.
#
#   "mouse"   — Mouse title-case symbols (Col1a1 -> COL1A1). Works for
#               >95% of mouse-human gene pairs: most symbols are identical
#               across species, differing only in capitalization. Genes with
#               genuine symbol divergence (some immune/MHC/cluster loci) will
#               fall out at the intersection step.
#               TODO: upgrade to full mouse-human ortholog mapping via
#               org.Mm.eg.db when mouse digit-tip dataset is integrated.
#
#   "xenopus" — Strip .L/.S (and optional .N) homeolog suffix, then
#               toupper(). Both subgenome copies of the same gene map to
#               one human symbol; use aggregate_to_hgnc() to sum/average
#               before calling this if homeologs haven't been combined yet.
#
# Dependencies:
#   AnnotationDbi (always loaded with Seurat ecosystem)
#   org.Hs.eg.db — for "human" mode alias lookup
#   Install: BiocManager::install(c("org.Hs.eg.db", "org.Mm.eg.db"))
#
# Functions:
#   harmonize_to_hgnc(gene_vec, species)
#   aggregate_to_hgnc(mat, species, filter_fn)
# =============================================================================

if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  stop("AnnotationDbi required: BiocManager::install('AnnotationDbi')")

# ---------------------------------------------------------------------------
# harmonize_to_hgnc()
#
# Maps a character vector of gene names to canonical HGNC uppercase symbols.
#
# Args:
#   gene_vec : character vector (e.g. rownames of an expression matrix)
#   species  : "human", "mouse", or "xenopus" (see file header)
#
# Returns:
#   Named character vector, length == length(gene_vec).
#   Names  = original gene names.
#   Values = HGNC uppercase symbol, or the uppercased input if no alias
#            match is found (never NA for "mouse"/"xenopus").
# ---------------------------------------------------------------------------
harmonize_to_hgnc <- function(gene_vec,
                               species = c("human", "mouse", "xenopus")) {
  species <- match.arg(species)

  # --- Xenopus: strip homeolog suffix ---
  if (species == "xenopus") {
    stripped <- toupper(sub("\\.(L|S)(\\.[0-9]+)?$", "", gene_vec))
    return(setNames(stripped, gene_vec))
  }

  # --- Mouse: uppercase (same symbol, different capitalization convention) ---
  if (species == "mouse") {
    return(setNames(toupper(gene_vec), gene_vec))
  }

  # --- Human / axolotl: org.Hs.eg.db ALIAS -> official HGNC SYMBOL ---
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("org.Hs.eg.db required: BiocManager::install('org.Hs.eg.db')")

  upper <- toupper(gene_vec)

  suppressMessages(
    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys      = unique(upper),
      column    = "SYMBOL",
      keytype   = "ALIAS",
      multiVals = "first"
    )
  )

  # mapped is a named vector: alias -> official SYMBOL (NA if not found)
  resolved <- unname(mapped[upper])

  # Fall back to the uppercased input for any name not in the alias table.
  # This preserves genes that are valid HGNC symbols but happen not to
  # appear as explicit aliases (very new genes, axolotl-specific names).
  result <- ifelse(is.na(resolved), upper, resolved)
  setNames(result, gene_vec)
}

# ---------------------------------------------------------------------------
# aggregate_to_hgnc()
#
# Harmonizes the rownames of a dense matrix and collapses duplicate HGNC
# symbols by averaging — needed when multiple input rows resolve to the same
# symbol (e.g., two SwissProt abbreviations for the same gene; Xenopus
# homeologs that weren't summed upstream).
#
# Args:
#   mat       : numeric matrix, rows = genes, cols = conditions/timepoints
#   species   : passed to harmonize_to_hgnc()
#   filter_fn : optional function(gene_vec) -> logical; TRUE = keep before
#               harmonizing (e.g. remove Trinity IDs before alias lookup)
#
# Returns:
#   Dense matrix, rows = unique HGNC symbols, cols = same as input.
# ---------------------------------------------------------------------------
aggregate_to_hgnc <- function(mat, species, filter_fn = NULL) {
  genes <- rownames(mat)

  if (!is.null(filter_fn)) {
    keep  <- filter_fn(genes)
    mat   <- mat[keep, , drop = FALSE]
    genes <- rownames(mat)
  }

  sym_map <- harmonize_to_hgnc(genes, species)

  # Drop any genes that resolved to empty string (shouldn't happen, but guard)
  valid <- nzchar(sym_map)
  mat   <- mat[valid, , drop = FALSE]
  syms  <- sym_map[valid]

  # rowsum() sums each column group; dividing by group size gives the mean.
  # Faster than a for-loop over unique symbols.
  agg_sum   <- rowsum(mat, group = syms, reorder = FALSE)
  grp_count <- tabulate(match(syms, rownames(agg_sum)))
  sweep(agg_sum, 1, grp_count, FUN = "/")
}
