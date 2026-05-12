# Narrative — Dataset Integration Basics

## Audience
Lab members who know the biology deeply — regeneration, BMS, axolotl anatomy — but are not bioinformaticians. They know what bulk RNA-seq is and can read a heatmap, but do not know what a UMAP is, what clustering algorithms do, why batch correction is needed, or what "pseudobulk" means. The goal is to build intuition for the bioinformatics steps in biological terms, not to teach them to run the code.

## One sentence
Single-cell RNA-seq integration is a sequence of specific technical fixes — each one solves a problem that would otherwise produce misleading results — and understanding each step helps you trust the output and ask better questions with it.

## Core questions
1. What does the raw data look like, and why can't we just use it directly?
2. What does each processing step actually do, and what would go wrong if we skipped it?
3. How do we get from three species' worth of single-cell data to a single comparison against our bulk RNA-seq?

## Key findings
- Each integration step solves one specific problem: QC removes noise, normalization removes depth bias, Harmony removes batch effects, annotation gives clusters biological meaning, ortholog mapping makes gene names compatible, pseudobulk converts sparse single-cell data into something bulk-comparable.
- The final result — a Spearman correlation heatmap — is only interpretable because each upstream step was done correctly. Skipping any step would change the answer.

## Emphasis
- **Foreground**: The flow — what each step does, what it fixes, what would happen without it. Use simple analogies throughout.
- **Foreground**: The "why this order matters" — each step produces input for the next; you cannot reorder them.
- **Background/one slide**: The actual biological result. This talk is about the method, not the finding.
- **Cut**: Per-dataset technical details, specific parameters, code. Those belong in the pipeline talk.
- **Cut**: Ortholog database comparisons, biomaRt vs toupper tradeoffs — too technical for this audience.

## Story arc

### Act 1: Why single-cell?
- Bulk RNA-seq gives you an average across millions of cells — it cannot tell you which cell types changed.
- Single-cell RNA-seq measures each cell individually — you can see which populations respond to a treatment.
- But scRNA-seq data is noisy, sparse, and technically complex — you can't just load it and start analyzing.

### Act 2: What the raw data looks like
- A count matrix: genes (rows) × cells (columns), mostly zeros (a typical cell only "sees" 10-30% of its genes in any one capture).
- Multiple samples per dataset, multiple species, different gene naming conventions.
- Before we can compare anything, we need to clean, normalize, and harmonize.

### Act 3: The six integration steps
One slide per step, each answering: what's the problem, what's the fix, what does the QC check look like.
1. QC — remove broken/fake cells
2. Normalization — put all cells on the same scale
3. Harmony — remove sequencing run noise
4. Cell type annotation — give clusters biological identity
5. Ortholog mapping — translate gene names across species
6. Pseudobulk — collapse single cells into a bulk-comparable profile

### Act 4: The comparison
- Spearman correlation: rank genes by expression, ask which scRNA-seq timepoint ranks them most similarly to our bulk data.
- The heatmap: one number per bulk condition × scRNA-seq timepoint pair.
- The answer: Xenopus froglet 7-14 dpa is the closest match.
