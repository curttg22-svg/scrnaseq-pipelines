# =============================================================================
# 00 — Download raw data from GEO (GSE165901)
# Run once from the xenopus/ project root: source("R/00_download_data.R")
#
# GSE165901 contains two experimental groups:
#   BL  (blastema) — froglet post-amputation, NON-regenerative
#   LBst (limb bud stage) — tadpole NF stages, REGENERATIVE window
#
# Excluded samples (rationale):
#   GSM5057666  Xen_Pool_LBst54_BL0dpa  — mixed regen + non-regen in one library
#   GSM5057667  Xen_Transplant           — tissue transplant experiment
#   GSM5045048  Axo_BL11dpa_rep1        — axolotl, not Xenopus
#   Axo_BL11dpa_rep2                    — axolotl rep2, not Xenopus
# =============================================================================

library(GEOquery)
library(tools)

raw_dir  <- file.path("data", "raw")
anno_dir <- file.path("data", "annotation")
dir.create(raw_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(anno_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Froglet BL (blastema) samples — non-regenerative, post-amputation
# ---------------------------------------------------------------------------

# Original samples (GSM5045xxx series, animal IDs vary)
SAMPLES_BL_ORIGINAL <- c(
  "GSM5045045" = "Xen_BL0dpa",
  "GSM5045046" = "Xen_BL3dpa",
  "GSM5045047" = "Xen_Pool_BL7_10_14dpa"
)

# Extended BL samples (GSM5057xxx series)
# GSM5057665 and GSM5057660 were added in the first extension pass.
# GSM5057655/656 are FACS-sorted from animal 107606 (same animal as the
# pool and late samples), supplementing the 0 dpa and 3 dpa coverage.
SAMPLES_BL_EXTENDED <- c(
  "GSM5057655" = "Xen_FACS_BL0dpa",
  "GSM5057656" = "Xen_FACS_BL3dpa",
  "GSM5057665" = "Xen_BL14dpa",
  "GSM5057660" = "Xen_Pool_BL14_20_52dpa"
)

# ---------------------------------------------------------------------------
# Tadpole LBst (limb bud stage) samples — regenerative window
# NF stages 50-52: hindlimb bud developing, amputation-competent for regen.
# ---------------------------------------------------------------------------
SAMPLES_LBST <- c(
  "GSM5057657" = "Xen_LBst50",
  "GSM5057658" = "Xen_LBst51",
  "GSM5057659" = "Xen_LBst52"
)

# ---------------------------------------------------------------------------
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

download_all <- function(sample_map, raw_dir) {
  for (acc in names(sample_map))
    download_sample(acc, sample_map[[acc]], raw_dir)
}

message("=== Downloading GSE165901 Xenopus samples ===")

message("\n--- Froglet BL (original) ---")
download_all(SAMPLES_BL_ORIGINAL, raw_dir)

message("\n--- Froglet BL (extended/FACS supplement) ---")
download_all(SAMPLES_BL_EXTENDED, raw_dir)

message("\n--- Tadpole LBst (regenerative) ---")
download_all(SAMPLES_LBST, raw_dir)

message("\n=== Download complete ===")
message("Run 01_xenopus_qc_clustering.R next.")
