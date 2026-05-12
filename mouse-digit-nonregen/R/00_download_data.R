# =============================================================================
# 00 — Download raw data from GEO (GSE135985)
# Run once from the mouse-digit-nonregen/ project root: source("R/00_download_data.R")
#
# Storer et al. 2020 (Dev Cell) — PMID 31902657
# "Acquisition of a unique mesenchymal precursor-like blastema state
#  underlies successful adult mammalian digit tip regeneration"
#
# Same lab (Storer/Miller) as GSE143888 but earlier paper. Crucially, this
# dataset contains both REGENERATIVE and NON-REGENERATIVE amputations in
# the same experiment, enabling a direct within-study comparison.
#
# Regenerative: distal (tip) amputation — digit tip regrows
# Non-regenerative: proximal amputation — wound heals by fibrosis, no regrowth
#
# Samples selected:
#   Regenerative:    control + 7, 10, 14, 28, 56 dpa (6 samples, 2 reps at 14dpa)
#   Non-regenerative: 10 and 14 dpa  (3 samples: 2 reps at 10dpa, 1 at 14dpa)
#
# Excluded: embryonic (E11, E14), neonatal (P3), transgenic Dmp1CreERT2 samples
# =============================================================================

library(GEOquery)

raw_dir <- file.path("data", "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

SAMPLES <- c(
  # --- Regenerative (distal amputation) ---
  "GSM4038976" = "Regen_7dpa",
  "GSM4038977" = "Regen_10dpa",
  "GSM4038978" = "Regen_14dpa_A",
  "GSM4038979" = "Regen_14dpa_B",
  "GSM4038980" = "Regen_28dpa",
  "GSM4038981" = "Regen_56dpa",
  "GSM4038982" = "Regen_uninjured_A",
  "GSM4038983" = "Regen_uninjured_B",

  # --- Non-regenerative (proximal amputation) ---
  "GSM4038984" = "NonRegen_10dpa_A",
  "GSM4038985" = "NonRegen_10dpa_B",
  "GSM4038986" = "NonRegen_14dpa"
)

download_sample <- function(acc, dest_name, raw_dir) {
  dest <- file.path(raw_dir, dest_name)
  if (all(file.exists(file.path(dest, c("barcodes.tsv.gz",
                                         "features.tsv.gz",
                                         "matrix.mtx.gz"))))) {
    message("  Already present, skipping: ", dest_name)
    return(invisible(NULL))
  }
  dir.create(dest, showWarnings = FALSE)
  message("\nDownloading ", acc, " -> ", dest_name)
  supp_files <- getGEOSuppFiles(acc, makeDirectory = FALSE, baseDir = dest)
  for (f in rownames(supp_files)) {
    bn <- basename(f)
    if (grepl("matrix",  bn, ignore.case = TRUE)) file.rename(f, file.path(dest, "matrix.mtx.gz"))
    if (grepl("barcode", bn, ignore.case = TRUE)) file.rename(f, file.path(dest, "barcodes.tsv.gz"))
    if (grepl("feature|gene", bn, ignore.case = TRUE)) file.rename(f, file.path(dest, "features.tsv.gz"))
  }
  message("  Done: ", paste(list.files(dest), collapse = ", "))
}

message("=== Downloading GSE135985 (Storer 2020 regenerative + non-regenerative) ===")
for (acc in names(SAMPLES))
  download_sample(acc, SAMPLES[[acc]], raw_dir)

message("\n=== Download complete ===")
message("Run 01_mouse_nonregen_clustering.R next.")
message("NOTE: Mouse 10x gene IDs are usually title-case symbols (Col1a1).")
