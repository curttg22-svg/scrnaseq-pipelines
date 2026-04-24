# =============================================================================
# 00 — Download raw data from GEO (GSE165901)
# Run once from the project root: source("R/00_download_data.R")
# =============================================================================

library(GEOquery)
library(tools)

# Paths relative to project root
raw_dir  <- file.path("data", "raw")
anno_dir <- file.path("data", "annotation")
dir.create(raw_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(anno_dir, showWarnings = FALSE, recursive = TRUE)

# Sample map: GEO accession -> local folder name
# Original samples (GSM5045xxx series)
SAMPLES <- c(
  "GSM5045045" = "Xen_BL0dpa",
  "GSM5045046" = "Xen_BL3dpa",
  "GSM5045047" = "Xen_Pool_BL7_10_14dpa",
  "GSM5045048" = "Axo_BL11dpa_rep1"
  # GSM5045049 (axolotl rep2) intentionally excluded — severe batch effect
)

# Extended Xenopus regeneration timepoints (GSM5057xxx series)
# Added to capture late regeneration window (14-52 dpa) where Hedgehog
# pathway GOI are expected to show activity if present in Xenopus.
# Axolotl samples from this series excluded — CT fate-map enriched (~97%
# mesenchymal), not suitable for cell-type-balanced cross-species comparison.
SAMPLES_EXTENDED <- c(
  "GSM5057665" = "Xen_BL14dpa",
  "GSM5057660" = "Xen_Pool_BL14_20_52dpa"
)

# Set to TRUE to also download the original samples (skip if already present)
DOWNLOAD_ORIGINAL <- FALSE

download_samples <- function(sample_map, raw_dir) {
  for (acc in names(sample_map)) {
    dest <- file.path(raw_dir, sample_map[acc])

    # Skip if already downloaded
    if (all(file.exists(file.path(dest, c("barcodes.tsv.gz",
                                           "features.tsv.gz",
                                           "matrix.mtx.gz"))))) {
      message("  Already present, skipping: ", sample_map[acc])
      next
    }

    dir.create(dest, showWarnings = FALSE)
    message("\nDownloading ", acc, " -> ", dest)

    supp_files <- getGEOSuppFiles(acc, makeDirectory = FALSE, baseDir = dest)

    dl_files <- rownames(supp_files)
    for (f in dl_files) {
      bname <- basename(f)
      if (grepl("matrix", bname, ignore.case = TRUE))
        file.rename(f, file.path(dest, "matrix.mtx.gz"))
      if (grepl("barcode", bname, ignore.case = TRUE))
        file.rename(f, file.path(dest, "barcodes.tsv.gz"))
      if (grepl("feature|gene", bname, ignore.case = TRUE))
        file.rename(f, file.path(dest, "features.tsv.gz"))
    }
    message("  Done: ", paste(list.files(dest), collapse = ", "))
  }
}

message("Fetching GSE165901 supplementary file list...")

if (DOWNLOAD_ORIGINAL) {
  message("\n--- Original samples ---")
  download_samples(SAMPLES, raw_dir)
}

message("\n--- Extended Xenopus regeneration timepoints ---")
download_samples(SAMPLES_EXTENDED, raw_dir)

message("\n=== Download complete ===")
message("Now download the axolotl GTF annotation manually:")
message("  https://www.axolotlomics.org/resources")
message("  Save as: data/annotation/AmexT_v47-AmexG_v6.0-DD.gtf")
