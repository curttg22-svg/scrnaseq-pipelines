# =============================================================================
# 04 — Axolotl GTF annotation: AMEX60DD ID -> human gene symbol bridge
# Input:  data/annotation/AmexT_v47-AmexG_v6.0-DD.gtf
# Output: results/final_bridge.rds
#
# The axolotlomics GTF encodes gene names in compound format:
#   "AXOLOTLNAME [NR]|HUMANNAME [HS]"  — extract the [HS] portion
#   "GENENAME"                          — use as-is (already clean)
#   "AMEX60DD000005"                    — no annotation, skip
# =============================================================================

library(stringr)
library(dplyr)

results_dir <- "results"
gtf_path    <- file.path("data","annotation","AmexT_v47-AmexG_v6.0-DD.gtf")

if (!file.exists(gtf_path)) {
  stop("GTF not found at: ", gtf_path,
       "\nDownload from axolotlomics.org and place at the path above.")
}

# =============================================================================
# 1. PARSE GTF — gene-level rows only for speed
# =============================================================================

message("Parsing GTF (gene rows only)...")
gtf_raw <- read.table(
  pipe(paste0("grep -v '^#' ", shQuote(gtf_path),
              " | awk -F'\\t' '$3==\"gene\"'")),
  sep = "\t", header = FALSE, quote = "", comment.char = "",
  col.names = c("seqname","source","feature","start","end",
                "score","strand","frame","attributes")
)
message("Gene rows in GTF: ", nrow(gtf_raw))

# Extract gene_id (AMEX ID) and gene_name
gtf_raw$amex_id   <- str_match(gtf_raw$attributes,
                                 'gene_id\\s+"([^"]+)"')[, 2]
gtf_raw$gene_name <- str_match(gtf_raw$attributes,
                                 'gene_name\\s+"([^"]+)"')[, 2]

# =============================================================================
# 2. PARSE COMPOUND GENE NAMES
# =============================================================================

gtf_map <- gtf_raw |>
  mutate(
    # Extract [HS] human symbol if present
    hs_symbol = str_match(gene_name, "([^|]+)\\s*\\[HS\\]")[, 2],
    gene_symbol_clean = case_when(
      !is.na(hs_symbol)                                   ~ trimws(hs_symbol),
      !grepl("\\[NR\\]|AMEX60DD|LOC[0-9]", gene_name)    ~ trimws(gene_name),
      TRUE                                                 ~ NA_character_
    ),
    gene_symbol_clean = toupper(gene_symbol_clean)
  ) |>
  filter(!is.na(gene_symbol_clean), !is.na(amex_id)) |>
  select(amex_id, gene_symbol = gene_symbol_clean) |>
  distinct(amex_id, .keep_all = TRUE)

message("AMEX IDs with clean gene symbol: ", nrow(gtf_map))
# Expected: ~24,811 (2 unmapped are ERCC spike-ins)

# =============================================================================
# 3. FINALISE BRIDGE TABLE
# Renames columns to the downstream convention (amex_id, amex_symbol).
# Cross-species ortholog extension (AMEX -> Xenopus homeolog) would require
# an ENSEMBL ortholog table and is not implemented here — human symbol is used
# as the common key in 06_goi_visualization.R instead.
# =============================================================================

final_bridge <- gtf_map |>
  rename(amex_id = amex_id, amex_symbol = gene_symbol)

saveRDS(final_bridge, file.path(results_dir, "final_bridge.rds"))
message("Saved results/final_bridge.rds")
message("Sample rows:")
print(head(final_bridge, 10))
