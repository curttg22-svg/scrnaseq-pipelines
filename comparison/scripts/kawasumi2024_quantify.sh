#!/usr/bin/env bash
# =============================================================================
# kawasumi2024_quantify.sh
#
# Quantifies Kawasumi-Kita et al. 2024 (PRJNA785721) bulk RNA-seq using Salmon.
# Processes 24 SRR accessions across 5 sample groups, one at a time with
# immediate FASTQ cleanup to stay within ~13 GB free disk.
#
# Target samples:
#   Froglet controls (PAIRED-end, 150 bp, ~1.8 GB/run):
#     Control_7dpa  : SRR27895473-476  (4 replicates)
#     Control_14dpa : SRR27908757-760  (4 replicates)
#     Control_21dpa : SRR27911609-612  (4 replicates)
#
#   Tadpole regen (SINGLE-end, 78 bp, ~350 MB/run):
#     Regen_5dpa    : SRR17108583-585,SRR17108587-589  (3 bio reps x 2 fractions)
#     Regen_7dpa    : SRR17108577-582                  (3 bio reps x 2 fractions)
#
# PRE-FLIGHT (run before starting):
#   (a) Inputs  : SRA reads (streamed via prefetch); Ensembl XenLae2 cDNA FASTA
#   (b) Disk    : ~130 MB FASTA + ~1 GB index + <=5 GB peak per run (freed after)
#                 Net kept on completion: ~1.5 GB index + ~240 MB quant outputs
#   (c) Time    : ~30 min setup; ~5-15 min per SRR x 24 = 2-6 hrs total
#   (d) Existing: None found (confirmed May 2026)
#   (e) Minimal : Process one SRR at a time, delete FASTQs after quant
#
# Run from comparison/ directory:
#   bash scripts/kawasumi2024_quantify.sh
#
# Requires: salmon (bioconda), sra-tools (prefetch + fasterq-dump), pigz or gzip
# Install:  conda install -c bioconda -c conda-forge salmon sra-tools
# =============================================================================

set -euo pipefail

# --- Paths -------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")/data/xenopus_kawasumi2024"
REF_DIR="$BASE_DIR/ref"
IDX_DIR="$BASE_DIR/salmon_index"
QUANT_DIR="$BASE_DIR/quant"
TMP_DIR="$BASE_DIR/tmp_fastq"
TX2GENE="$REF_DIR/tx2gene.tsv"
TRANSCRIPTOME_URL="https://ftp.ensembl.org/pub/release-112/fasta/xenopus_laevis/cdna/Xenopus_laevis.XenLae2.cdna.all.fa.gz"
TRANSCRIPTOME_FA="$REF_DIR/Xenopus_laevis.XenLae2.cdna.all.fa.gz"

mkdir -p "$REF_DIR" "$IDX_DIR" "$QUANT_DIR" "$TMP_DIR"

# --- Sample manifest ---------------------------------------------------------
# Format: SRR_ID|sample_group|library_layout
declare -a SAMPLES=(
  "SRR27895473|Control_7dpa|PAIRED"
  "SRR27895474|Control_7dpa|PAIRED"
  "SRR27895475|Control_7dpa|PAIRED"
  "SRR27895476|Control_7dpa|PAIRED"
  "SRR27908757|Control_14dpa|PAIRED"
  "SRR27908758|Control_14dpa|PAIRED"
  "SRR27908759|Control_14dpa|PAIRED"
  "SRR27908760|Control_14dpa|PAIRED"
  "SRR27911609|Control_21dpa|PAIRED"
  "SRR27911610|Control_21dpa|PAIRED"
  "SRR27911611|Control_21dpa|PAIRED"
  "SRR27911612|Control_21dpa|PAIRED"
  "SRR17108583|Regen_5dpa|SINGLE"
  "SRR17108584|Regen_5dpa|SINGLE"
  "SRR17108585|Regen_5dpa|SINGLE"
  "SRR17108587|Regen_5dpa|SINGLE"
  "SRR17108588|Regen_5dpa|SINGLE"
  "SRR17108589|Regen_5dpa|SINGLE"
  "SRR17108577|Regen_7dpa|SINGLE"
  "SRR17108578|Regen_7dpa|SINGLE"
  "SRR17108579|Regen_7dpa|SINGLE"
  "SRR17108580|Regen_7dpa|SINGLE"
  "SRR17108581|Regen_7dpa|SINGLE"
  "SRR17108582|Regen_7dpa|SINGLE"
)

# --- Tool checks -------------------------------------------------------------
for tool in salmon prefetch fasterq-dump; do
  if ! command -v "$tool" &>/dev/null; then
    echo "ERROR: '$tool' not found."
    echo "Install with: conda install -c bioconda -c conda-forge salmon sra-tools"
    exit 1
  fi
done
echo "Tools OK: salmon $(salmon --version 2>&1 | head -1), prefetch $(prefetch --version 2>&1 | head -1)"

# --- Step 1: Download transcriptome ------------------------------------------
if [[ ! -f "$TRANSCRIPTOME_FA" ]]; then
  echo "Downloading X. laevis cDNA FASTA (~130 MB)..."
  curl -L "$TRANSCRIPTOME_URL" -o "$TRANSCRIPTOME_FA"
else
  echo "Transcriptome already present: $TRANSCRIPTOME_FA"
fi

# --- Step 2: Build tx2gene table from FASTA headers --------------------------
if [[ ! -f "$TX2GENE" ]]; then
  echo "Building tx2gene table from FASTA headers..."
  zcat "$TRANSCRIPTOME_FA" | awk '
    /^>/ {
      tx = ""; gid = ""; sym = ""
      # transcript ID is the first field after >
      match($0, />([^ ]+)/, a);  tx = a[1]
      # gene ID
      match($0, /gene:([^ ]+)/, b); gid = b[1]
      # gene symbol (may be absent for some transcripts)
      if (match($0, /gene_symbol:([^ ]+)/, c)) {
        sym = c[1]
      } else {
        sym = gid
      }
      print tx "\t" gid "\t" sym
    }
  ' > "$TX2GENE"
  echo "tx2gene: $(wc -l < "$TX2GENE") transcripts"
else
  echo "tx2gene already present: $TX2GENE"
fi

# --- Step 3: Build Salmon index ----------------------------------------------
if [[ ! -d "$IDX_DIR/info.json" ]] && [[ ! -f "$IDX_DIR/info.json" ]]; then
  echo "Building Salmon index (this takes ~5-10 min)..."
  salmon index \
    --transcripts "$TRANSCRIPTOME_FA" \
    --index "$IDX_DIR" \
    --threads 4 \
    --gencode
  echo "Index built: $IDX_DIR"
else
  echo "Salmon index already present: $IDX_DIR"
fi

# --- Step 4: Quantify each SRR -----------------------------------------------
N_TOTAL=${#SAMPLES[@]}
N_DONE=0

for entry in "${SAMPLES[@]}"; do
  IFS="|" read -r SRR GROUP LAYOUT <<< "$entry"
  OUT="$QUANT_DIR/$SRR"

  if [[ -f "$OUT/quant.sf" ]]; then
    echo "[$((++N_DONE))/$N_TOTAL] $SRR already quantified, skipping."
    continue
  fi

  echo ""
  echo "[$((N_DONE+1))/$N_TOTAL] Processing $SRR ($GROUP, $LAYOUT)..."

  # Download .sra file
  echo "  Prefetching $SRR..."
  prefetch "$SRR" --output-directory "$TMP_DIR" --max-size 10G

  # Convert to FASTQ
  echo "  Extracting FASTQ..."
  fasterq-dump \
    --split-files \
    --outdir "$TMP_DIR" \
    --temp "$TMP_DIR" \
    --threads 4 \
    "$TMP_DIR/$SRR/$SRR.sra"

  # Quantify
  echo "  Running Salmon quant..."
  if [[ "$LAYOUT" == "PAIRED" ]]; then
    salmon quant \
      --index "$IDX_DIR" \
      --libType A \
      --mates1 "$TMP_DIR/${SRR}_1.fastq" \
      --mates2 "$TMP_DIR/${SRR}_2.fastq" \
      --threads 4 \
      --validateMappings \
      --output "$OUT"
  else
    salmon quant \
      --index "$IDX_DIR" \
      --libType A \
      --unmatedReads "$TMP_DIR/${SRR}.fastq" \
      --threads 4 \
      --validateMappings \
      --output "$OUT"
  fi

  # Clean up immediately to free disk
  echo "  Cleaning up FASTQs..."
  rm -rf "$TMP_DIR/$SRR" "$TMP_DIR/${SRR}.fastq" \
         "$TMP_DIR/${SRR}_1.fastq" "$TMP_DIR/${SRR}_2.fastq"

  N_DONE=$((N_DONE + 1))
  echo "  Done. Disk free: $(df -h . | awk 'NR==2{print $4}')"
done

echo ""
echo "All $N_DONE/$N_TOTAL SRRs quantified."
echo "quant.sf files in: $QUANT_DIR/"
echo "Next step: source('R/project_kawasumi_mds.R') in R"
