# =============================================================================
# 00 — Download raw data from GEO (GSE143888)
# Run once from the mouse-digit/ project root: source("R/00_download_data.R")
#
# Johnson, Masias & Lehoczky 2020 (Dev Cell) — PMID 32097654
# "Cellular Heterogeneity and Lineage Restriction during Mouse Digit Tip
#  Regeneration at Single-Cell Resolution"
# Lehoczky lab, Brigham and Women's Hospital / Harvard
# Adult mouse digit tips, distal phalanx amputation model
# Regenerative: mouse digit tip CAN regenerate at the distal phalanx
#
# 5 samples: unamputated control + 4 post-amputation timepoints
# =============================================================================

library(GEOquery)

raw_dir <- file.path("data", "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

SAMPLES <- c(
  "GSM4276219" = "Mouse_0dpa",   # unamputated control
  "GSM4276220" = "Mouse_11dpa",  # 11 days post-amputation
  "GSM4276221" = "Mouse_12dpa",  # 12 days post-amputation
  "GSM4276222" = "Mouse_14dpa",  # 14 days post-amputation
  "GSM4276223" = "Mouse_17dpa"   # 17 days post-amputation
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

message("=== Downloading GSE143888 Mouse digit tip samples ===")
for (acc in names(SAMPLES))
  download_sample(acc, SAMPLES[[acc]], raw_dir)

message("\n=== Download complete ===")
message("Run 01_mouse_qc_clustering.R next.")
message("NOTE: Check gene ID format in features.tsv.gz after download.")
message("      Mouse 10x data uses title-case symbols (Col1a1) or Ensembl IDs.")
