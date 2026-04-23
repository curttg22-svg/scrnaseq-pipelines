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
SAMPLES <- c(
  "GSM5045045" = "Xen_BL0dpa",
  "GSM5045046" = "Xen_BL3dpa",
  "GSM5045047" = "Xen_Pool_BL7_10_14dpa",
  "GSM5045048" = "Axo_BL11dpa_rep1"
  # GSM5045049 (axolotl rep2) intentionally excluded — severe batch effect
)

message("Fetching GSE165901 supplementary file list...")
gse <- getGEO("GSE165901", GSEMatrix = FALSE)

for (acc in names(SAMPLES)) {
  dest <- file.path(raw_dir, SAMPLES[acc])
  dir.create(dest, showWarnings = FALSE)

  message("\nDownloading ", acc, " -> ", dest)
  gsm <- getGEO(acc)

  # GEO supplementary files for 10x MTX samples
  supp_files <- getGEOSuppFiles(acc, makeDirectory = FALSE, baseDir = dest)

  # Rename to standard 10x names if needed
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

message("\n=== Download complete ===")
message("Now download the axolotl GTF annotation manually:")
message("  https://www.axolotlomics.org/resources")
message("  Save as: data/annotation/AmexT_v47-AmexG_v6.0-DD.gtf")
