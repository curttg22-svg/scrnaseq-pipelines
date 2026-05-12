# Presentation Slide Outline

## Metadata
- **Title**: Single-cell RNA-seq integration — the basics
- **Date**: May 2026
- **Author**: T. Curtis <curttg22@wfu.edu>
- **Source**: `comparison/bulk_vs_scrna_correlation.R`, `xenopus/R/`, `mouse-digit/R/`
- **Presentation Skill Version**: 1.0

---

## Slides

---

### Slide 1: Title slide
**Title**: Single-cell RNA-seq integration\nHow the data goes from raw reads to a cross-species comparison

**Type**: title

**Content**:
- Lab meeting | May 2026

---

### Slide 2: Bulk RNA-seq gives you a population average — single-cell shows you each cell individually
**Title**: Bulk RNA-seq averages across all cells; single-cell RNA-seq measures each one separately

**Content**:
- **Bulk RNA-seq**: homogenize tissue, extract RNA from millions of cells together — you get the average gene expression of the whole population, but you cannot see which cell types are responsible for any given change
- **Single-cell RNA-seq**: dissociate tissue into a cell suspension, capture each cell individually in a microfluidic droplet, sequence the RNA from that one cell — you get a gene expression profile per cell, and can then group cells by type
- The tradeoff: scRNA-seq is noisier (each cell is sampled much more shallowly than bulk), and the data is far more complex to process

**Notes**: Use the analogy of a smoothie vs individual fruits. Bulk RNA-seq is like blending everything and asking "what flavor is this?" Single-cell is like sorting each piece of fruit individually before blending. You learn much more about composition, but the measurement of each individual piece is less precise.

---

### Slide 3: A single-cell experiment produces a sparse matrix — most values are zero
**Title**: Raw single-cell data is a genes × cells matrix — mostly zeros due to shallow sequencing depth

**Content**:
- Each cell is sequenced shallowly — a typical cell captures expression from only 10-30% of its genes; the rest read as zero (not absent, just not detected)
- This sparsity is not a failure; it is a fundamental property of the technology — it means you cannot directly compare individual cells the way you compare bulk samples
- The matrix for one dataset: tens of thousands of genes (rows) × tens of thousands of cells (columns), with the vast majority of entries being zero

**Notes**: The sparsity point is important for building intuition for later steps. When the audience hears "pseudobulk" later, they need to understand why you can't compare a single cell to bulk data — most genes are zero in any one cell, so the correlation would be meaningless. Averaging across many cells restores enough signal to be bulk-comparable.

---

### Slide 4: Our three datasets span the regeneration spectrum across species
**Title**: Three species were chosen to cover the full range from regenerative to non-regenerative

**Content**:
- **Axolotl** (collaborator) — 43,083 cells — the gold standard of limb regeneration; forms a full blastema and regrows the limb completely
- **Xenopus** (Lin et al. 2021) — 65,911 cells — stage-dependent: tadpole regenerates limbs, adult froglet does not; gives us both outcomes from one paper
- **Mouse digit tip** (Johnson et al. 2020) — 37,743 cells — limited but real regeneration of the very tip; regeneration-competent under specific conditions

**Notes**: The species selection is not arbitrary — we want a reference that covers the full spectrum so that when our BMS bulk data matches one state over another, the result has biological meaning. If all three species were fully regenerative, a high correlation with any of them would be uninterpretable.

---

### Slide 5: Section divider — Step 1: Quality control
**Title**: Step 1: Quality control

**Type**: section

---

### Slide 6: 10x droplet capture is imperfect — QC identifies and removes three classes of bad data
**Title**: QC step: not every "cell" in the data is actually a cell — three categories of noise must be removed

**Figure**: [Xenopus pre-filter QC violin](figures/xenopus_qc_violin.png)

**Figure Position**: right

**Figure Caption**: Pre-filter QC violin. Thresholds separate real cells from noise.

**Content**:
- **Empty droplets**: ambient RNA, no cell — too few genes detected; removed by minimum gene threshold
- **Doublets**: two cells in one droplet — too many genes; removed by maximum gene threshold
- **Dying cells**: cytoplasm leaked, only mitochondrial RNA remains — flagged by high mito %

**Notes**: Analogy: running a gel and discarding bands in the wrong molecular weight range. The QC step is like the biologist's instinct to throw out a sample that looks wrong before putting it on the machine — here we are doing that computationally, after sequencing, using the data itself as the quality indicator.

---

### Slide 7: Section divider — Step 2: Normalization
**Title**: Step 2: Normalization

**Type**: section

---

### Slide 8: Cells are sequenced to different depths — normalization removes this technical bias
**Title**: Normalization puts all cells on the same scale so differences reflect biology, not sequencing depth

**Content**:
- Different cells end up with different total read counts depending on how much RNA they contained and how efficiently the library was prepared — a cell with 5,000 total UMI reads is not "more expressed" than one with 2,000; it was just captured more efficiently
- **Log-normalization**: divide each cell's gene counts by its total counts (scaling to a common depth), then apply a log transformation — after this, a value of 2 means the same thing in every cell regardless of how deeply it was sequenced
- Without normalization, cells with more total reads would appear to express everything more highly, and clustering would separate cells by sequencing depth rather than cell type

**Notes**: Analogy: if you counted words in documents of different lengths, a 10-page document would have more of every word than a 1-page document — even if they're about the same topic. Normalization divides by document length first so you can compare the relative frequency of words, not the absolute counts. The log transformation compresses the dynamic range (a gene expressed 100x more isn't 100x more biologically different — the log scale brings this into a comparable range).

---

### Slide 9: Section divider — Step 3: Dimensionality reduction and batch correction
**Title**: Step 3: UMAP and Harmony — visualize structure, remove batch effects

**Type**: section

---

### Slide 10: UMAP reduces 30,000 genes to two dimensions — nearby cells have similar transcriptomes
**Title**: UMAP is a map of transcriptional similarity — cells near each other express similar genes

**Figure**: [Xenopus UMAP annotated by cell type](figures/xenopus_umap_annotated.png)

**Figure Position**: right

**Figure Caption**: Xenopus UMAP. Each dot is one cell; proximity = transcriptional similarity.

**Content**:
- Each cell is described by ~30,000 gene values — UMAP compresses this to 2 dimensions
- Cells with similar transcriptomes land near each other; clusters of dots = groups with shared expression programs
- These clusters are the raw material for cell type annotation — the algorithm finds the groups; biology names them

**Notes**: UMAP (Uniform Manifold Approximation and Projection) is a nonlinear dimensionality reduction. It is NOT a coordinate system — distances between distant clusters are not meaningful, only local neighborhoods. The analogy: it's like a map of a city drawn by asking "who lives near whom?" — neighborhoods are correct but the scale across the whole map is distorted. The key thing the audience should take away is that the clusters are real (cells that are similar cluster together) even if the exact position on the plot is not physically meaningful.

---

### Slide 11: Harmony removes technical batch effects without erasing biological differences
**Title**: Harmony: 10x sequencing runs introduce noise — Harmony removes it while keeping biology intact

**Figure**: [Xenopus UMAP by timepoint](figures/xenopus_umap_timepoint.png)

**Figure Position**: right

**Figure Caption**: Post-Harmony UMAP by timepoint. Colors mixing within clusters = batch corrected.

**Content**:
- Each sequencing run introduces technical noise — without correction, cells cluster by run, not biology
- A fibroblast from timepoint 1 would cluster separately from a timepoint 2 fibroblast despite being the same cell type
- Harmony removes run-level noise; QC check: every cluster must contain cells from multiple timepoints

**Notes**: Analogy: if three different people pipette the same samples on three different days, there will be technical differences between their work even if the biology is identical. Batch correction is the computational equivalent of accounting for who ran which sample. The key intuition is that Harmony does NOT erase biological differences between timepoints — it specifically removes the variance that correlates with "which run was this" while preserving the variance that correlates with "what time point was this."

---

### Slide 12: Section divider — Step 4: Cell type annotation
**Title**: Step 4: Cell type annotation

**Type**: section

---

### Slide 13: Clustering groups similar cells — annotation tells you what those groups actually are
**Title**: Clustering finds groups of transcriptionally similar cells — annotation assigns biological identity to each group

**Figure**: [Xenopus canonical marker dotplot](figures/xenopus_marker_dotplot.png)

**Figure Position**: right

**Figure Caption**: Rows = clusters; columns = marker genes. Dot size = % expressing; color = level.

**Content**:
- Clustering is purely mathematical — it produces numbered groups with no biological identity yet
- Annotation: if a cluster expresses col1a1 + vim but not epcam or ptprc, it is likely a fibroblast
- The dotplot is the audit trail — documents why each label was assigned

**Notes**: This step requires biological expertise, not just computation. An algorithm can identify that Cluster 3 is transcriptionally distinct; it cannot tell you that col1a1+ vim+ cells are fibroblasts without knowledge of cell biology. This is why automated annotation tools are less reliable for non-standard species — there is no pre-built reference for axolotl limb cell types. A human expert who knows the marker genes must make the assignment.

---

### Slide 14: Section divider — Step 5: Cross-species ortholog mapping
**Title**: Step 5: Cross-species gene name harmonization

**Type**: section

---

### Slide 15: Each species uses different gene naming conventions — harmonization is required before any comparison
**Title**: "Ptch1", "ptch1.L", and "PTCH1" are the same gene written three different ways

**Content**:
- **Mouse**: title-case symbols — "Ptch1", "Col1a1", "Actb" — convention for mouse gene names
- **Xenopus laevis**: allotetraploid (genome duplicated twice) — every gene exists as two homeolog copies: "ptch1.L" (long chromosome) and "ptch1.S" (short chromosome) — both must be summed before comparison
- **Human/Axolotl**: HGNC uppercase — "PTCH1", "COL1A1", "ACTB" — standard for human gene databases; the axolotl collaborator dataset already uses this convention
- Without conversion, a string match would find zero shared genes between species — the pipeline converts all to HGNC uppercase before any comparison

**Notes**: The Xenopus homeolog situation is worth a sentence of explanation — Xenopus laevis underwent two rounds of whole-genome duplication in its evolutionary history, so it has two slightly diverged copies of almost every gene. These are called homeologs (.L and .S). When comparing to human or mouse (which have one copy), we sum the expression of both homeolog copies together to get the equivalent single-gene value. This is analogous to summing the two alleles in a diploid — we want total gene dosage, not per-allele dosage.

---

### Slide 16: Section divider — Step 6: Pseudobulk construction
**Title**: Step 6: Pseudobulk construction

**Type**: section

---

### Slide 17: Individual cells are too sparse to compare to bulk — averaging across cells restores the signal
**Title**: Pseudobulk: averaging all cells per timepoint converts sparse single-cell data into a bulk-comparable profile

**Content**:
- Recall that each cell detects only 10-30% of its genes — trying to correlate a single cell to a bulk RNA-seq profile would be comparing 2,000 genes to 30,000 genes; most would be zero in the cell by chance
- Pseudobulk: take all cells from a given timepoint and compute the mean expression across all of them — zero values average out, and the resulting profile is dense and bulk-comparable
- Result: a genes × timepoints matrix for each species, in the same format as our bulk RNA-seq data — now correlation is possible

**Notes**: The pseudobulk step is conceptually the inverse of why we use single-cell in the first place. We use single-cell to get cell type resolution, but then we collapse it for the bulk comparison. The reason is that our BMS bulk data represents a whole-tissue response — comparing it to pseudobulk (which also represents a whole-tissue average) is the right comparison. Comparing to one cell type's pseudobulk would answer a more specific question, which is a downstream analysis.

---

### Slide 18: Section divider — Step 7: Spearman correlation
**Title**: Step 7: Spearman rank correlation

**Type**: section

---

### Slide 19: Spearman rank correlation asks which scRNA-seq state has the most similar expression pattern to our bulk data
**Title**: Spearman correlation: rank genes by expression, then ask which scRNA-seq timepoint has the most similar ranking

**Content**:
- **Why not just compare numbers directly**: expression scales differ across species and assay types — a "high" value in mouse RNA-seq is not the same number as a "high" value in Xenopus scRNA-seq
- **Rank-based solution**: instead of comparing raw values, replace each gene's expression with its rank (1st most expressed, 2nd most expressed, etc.) within that profile — now every profile is on the same 1-to-N scale regardless of species or assay
- **The output**: one Spearman rho value per pair of bulk condition and scRNA-seq timepoint — rho near +1 means same genes are highly ranked in both, rho near 0 means no relationship

**Notes**: Analogy: instead of comparing test scores (which are on different scales between different exams), rank students by their score in each exam. Then ask: is the ranking of students similar across exams? Spearman does the same thing with genes — it asks whether the same genes tend to be relatively high-expressed (or low-expressed) in both the bulk and the scRNA-seq profile, regardless of the absolute values.

---

### Slide 20: The heatmap result: all bulk conditions most closely match Xenopus froglet 7-14 dpa
**Title**: The correlation heatmap reads like a table: warmer tiles mean more similar transcriptional programs

**Type**: figure

**Figure**: [Bulk vs scRNA-seq correlation heatmap](figures/heatmap_perpairwise.png)

**Figure Caption**: Spearman rho — bulk RNA-seq (rows) vs scRNA-seq pseudobulk (columns). Warmer = more similar.

**Notes**: Walk the audience through reading the heatmap before pointing to the result. Rows = our 4 bulk conditions (BMS top block, EtOH bottom). Columns = all scRNA-seq timepoints, grouped by species (Axolotl | Xenopus | Mouse). Each tile is one rho value. The warm patch at LinBL_7-14dpa is the signal. Xenopus froglet 7-14 dpa is a non-regenerative wound healing response — the digit closes but does not regrow. This is the most biologically similar state to our BMS bulk data.

---

### Slide 21: Each step is required — skipping any one changes the answer
**Title**: Each step solves one problem — skipping any step would change the result

**Content**:
- Skip QC → phantom clusters from empty droplets; annotation is unreliable
- Skip normalization → high-depth cells appear to "express more" everything; clusters separate by depth
- Skip Harmony → cells cluster by sequencing run, not cell type; annotation would be wrong
- Skip pseudobulk → individual cells too sparse to correlate with bulk; all rho near zero

**Notes**: This slide is the conceptual summary of why the pipeline is structured the way it is. Each step has a specific failure mode if skipped — and the failure modes compound (bad QC makes normalization worse; bad normalization makes batch correction harder). The order matters: you cannot correct batch effects before normalizing, and you cannot annotate cell types before removing batch effects.

---

### Slide 22: From raw reads to cross-species correlation — the full flow
**Title**: The full pipeline: six steps convert raw data into a cross-species comparison

**Content**:
- **QC → Normalization**: remove bad droplets; put all cells on a common scale
- **UMAP → Harmony**: visualize transcriptional similarity; remove sequencing batch effects
- **Annotation**: label clusters using canonical marker genes — biology names what the algorithm found
- **Ortholog mapping → Pseudobulk → Spearman**: align gene names, average cells per timepoint, correlate with bulk

**Notes**: This is the summary slide. The audience should leave with a mental model of why each step exists and what would go wrong without it. The result — BMS bulk data most closely resembles Xenopus froglet 7-14 dpa wound healing — is only interpretable because each upstream step was done correctly.

---
