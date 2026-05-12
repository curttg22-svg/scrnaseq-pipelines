# Narrative

## Audience
Lab members. They know regeneration biology and BMS experimental conditions. They are not deeply familiar with single-cell RNA-seq methods or why specific QC steps matter.

## One sentence
This pipeline systematically integrates multi-species regeneration scRNA-seq datasets into a validated compendium, then maps the lab's bulk RNA-seq to it via pseudobulk Spearman correlation — every step has a QC checkpoint that makes the result trustworthy.

## Core questions
1. How does the pipeline turn raw species-specific scRNA-seq datasets into something we can compare to our bulk data?
2. What checks ensure each step is accurate, not an artifact?
3. What does the initial result tell us about BMS biology?

## Key findings
- The pipeline is reproducible and validated at each step: QC filtering, batch correction, ortholog harmonization, and correlation strategy all have explicit checks.
- BMS conditions (Reamputation most strongly) consistently match Xenopus froglet 7-14 dpa (rho = 0.640) above any regenerative blastema — suggesting a non-regenerative wound healing transcriptomic signature.

## Emphasis
- **Foreground**: The pipeline architecture and QC logic — this is a lab resource talk, not a biology results talk.
- **Foreground**: The cross-species ortholog harmonization challenge and how we handle it (Trinity IDs, homeologs, common gene set).
- **Background/one slide**: The initial correlation finding as validation that the pipeline produces interpretable biology.
- **Cut/appendix**: Hedgehog pathway figures, individual tadpole comparisons, common gene set heatmap — these are follow-up analyses, not pipeline validation.

## Story arc

### Act 1: Motivation
- The lab has BMS bulk RNA-seq (4 conditions). We need to know what regenerative context this transcriptome resembles.
- Answer requires a curated, validated cross-species scRNA-seq compendium.

### Act 2: Per-dataset deep dive (NEW — expanded May 2026)
- **Axolotl collaborator** (private/unpublished): ShinyCell HDF5 format, pre-processed by collaborator, HGNC symbols, CCA UMAP pre-computed. Show cell type UMAP (15 types), condition UMAP (Regen/NonRegen separation), housekeeping QC violin (ACTG1+GAPDH stable across conditions).
- **Xenopus Lin 2021**: CellRanger → Seurat v5 → Harmony. Show pre-filter QC violin, timepoint UMAP (Harmony QC — interleaving), condition UMAP (BL vs LBst separation), annotated UMAP, canonical marker dotplot.
- **Mouse digit tip Johnson 2020**: Same pipeline. Show pre-filter QC violin, timepoint UMAP (Harmony QC), cluster UMAP (not yet annotated).

### Act 3: Cross-species harmonization
- Three strategies: direct HGNC (axolotl collab), suffix-strip toupper (Xenopus), toupper-only (mouse)
- Common gene set (4,545 genes) removes gene-count bias from rho comparison
- Per-pairwise vs common set strategies produce consistent rankings — validates robustness

### Act 4: Pseudobulk + correlation (condensed)
- Pseudobulk construction (mean log-norm per timepoint across all cells)
- Spearman rank correlation (why rank: scale-invariant across species)
- QC: two gene-set strategies — per-pairwise and common 4,545-gene set

### Act 3: Initial result
- The heatmap: all 4 bulk conditions correlate most strongly with Xenopus froglet 7-14 dpa.
- Reamputation zoom: best rho = 0.640 (Xenopus froglet 7-14dpa), axolotl and mouse digit lower.
- What this tells us about BMS, and what comes next.
