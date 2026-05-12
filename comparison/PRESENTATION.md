# Presentation Slide Outline

## Metadata
- **Title**: Multi-species scRNA-seq compendium pipeline
- **Date**: May 2026
- **Author**: T. Curtis <curttg22@wfu.edu>
- **Source**: `comparison/bulk_vs_scrna_correlation.R`, `xenopus/R/`, `mouse-digit/R/`, `comparison/R/`
- **Presentation Skill Version**: 1.0

---

## Slides

---

### Slide 1: Title slide
**Title**: A cross-species scRNA-seq compendium for contextualizing bulk RNA-seq\nPipeline walkthrough, QC validation, and conceptual rationale

**Type**: title

**Content**:
- Lab meeting | May 2026

---

### Slide 2: The experimental question — what transcriptional state does BMS induce?
**Title**: BMS changes gene expression, but we need a reference to interpret what state it creates

**Content**:
- The lab's BMS bulk RNA-seq (Acute, 24hrs, Pretreatment, Reamputation) measures a transcriptional response to NF-kB inhibition during axolotl limb amputation
- Bulk RNA-seq tells you what changed — it does not tell you which cell types drive those changes, or whether the resulting state resembles regeneration, wound healing, or something else
- Single-cell RNA-seq from species with known regenerative outcomes provides that reference: we can ask which published biological state our bulk data most resembles

**Notes**: This is the core conceptual motivation. Bulk RNA-seq is powerful but context-free — you see that 300 genes change, but without a reference frame you cannot say whether that looks like a regenerating blastema or a non-regenerating scar. The compendium provides that frame.

---

### Slide 3: The compendium spans three species at different positions on the regeneration spectrum
**Title**: Three species, four datasets, 145,000 cells spanning the regeneration spectrum

**Content**:
- **Axolotl** (collaborator, unpublished) — 43,083 cells | regen dpa3/14/23 + non-regen limb | 15 cell types | HGNC symbols, ShinyCell HDF5
- **Xenopus** (Lin et al. 2021) — 65,911 cells | froglet BL (non-regen, 5 dpa) + tadpole LBst (regen, NF50-52) | 21 cell types | 10x CellRanger, Seurat v5
- **Mouse digit tip** (Johnson et al. 2020, GSE143888) — 37,743 cells | 5 dpa timepoints | 10x CellRanger, Seurat v5

**Notes**: The species were chosen to span the regeneration spectrum: axolotl = full limb regeneration; Xenopus = stage-dependent (tadpole yes, froglet no); mouse digit tip = partial but real regrowth. This range means if BMS bulk data correlates with one state over others, the result has biological meaning — we are not just comparing to random transcriptomes.

---

### Slide 4: Section divider — Dataset 1: Collaborator Axolotl
**Title**: Dataset 1: Collaborator Axolotl

**Type**: section

---

### Slide 5: Collaborator data is pre-processed — we read it directly without re-clustering
**Title**: Using the collaborator's processed ShinyCell export avoids re-doing weeks of annotation work

**Content**:
- **Why**: Re-clustering and annotating 43,000 cells would require alignment to the axolotl genome, QC decisions, and weeks of marker-based annotation — the collaborator has already done this carefully for their unpublished study
- **Format**: ShinyCell HDF5 export — `merged_objassay_RNA.h5` (284 MB, dense genes × cells matrix), `merged_objgene.rds` (gene name → row index), `merged_objmeta.rds` (cell metadata with 15 annotated cell types)
- **Ortholog advantage**: genes are already HGNC symbols — no Trinity ID mapping, no suffix stripping, direct intersection with our bulk gene names

**Notes**: ShinyCell is a Shiny app export format from the Seurat ecosystem. It stores the normalized expression matrix as a dense HDF5 rather than a sparse matrix. We read it in 1,000-gene chunks to stay within RAM limits. The key benefit over getting an RDS file is that we don't need to load the full Seurat object — we read only the rows (genes) we need.

---

### Slide 6: The collaborator UMAP shows 15 annotated cell types across both conditions
**Title**: Collaborator axolotl: 15 annotated cell types, cleanly separated in the pre-computed UMAP

**Figure**: [Axolotl collaborator UMAP — cell types](figures/axo_collab_umap_celltypes.png)

**Figure Position**: right

**Figure Caption**: CCA UMAP, 43,083 cells, 15 annotated cell types (collaborator annotation).

**Content**:
- 43,083 cells: Regen 34,535 | NonRegen 12,548
- Regenerative Blastema Fibroblast: 3,832 Regen vs 24 NonRegen — regen-exclusive population
- Collaborator annotation used as-is — expert knowledge from the lab that generated the data

**Notes**: The CCA UMAP was computed by the collaborator. We trust it because it's from the group that generated the data and has biological expertise in this specific system. If we re-clustered from scratch without their guidance we would likely make annotation errors.

---

### Slide 7: Regen and non-regen cells are transcriptionally distinct — biological signal is preserved
**Title**: Regen and non-regen conditions separate in UMAP, confirming the dataset captures biology, not batch

**Figure**: [Axolotl collaborator UMAP — condition](figures/axo_collab_umap_condition.png)

**Figure Position**: right

**Content**:
- Blue = Regen (dpa3/14/23): 34,535 cells | Red = NonRegen (proximal limb): 12,548 cells
- Condition separation confirms the two states have distinct transcriptional programs
- Complete mixing would mean the dataset cannot distinguish regen from non-regen — it would be useless as a reference

**Notes**: This is a key QC check: condition separation should be biological, not batch. The Regenerative Blastema Fibroblast cluster is nearly absent in NonRegen (24 vs 3,832 cells) — this is not a batch artifact, it reflects the fact that blastema fibroblasts are a regeneration-specific cell state.

---

### Slide 8: Housekeeping gene stability confirms normalization is consistent across conditions
**Title**: ACTG1 and GAPDH are stable across all conditions — normalization artifacts would show as shifts

**Figure**: [Axolotl housekeeping QC violin](figures/axo_collab_housekeeping_violin.png)

**Type**: figure

**Figure Caption**: ACTG1 (actin proxy; ACTB absent) and GAPDH across Regen dpa3/14/23 and NonRegen. Log-norm.

**Notes**: Housekeeping gene QC is conceptually important: if log-normalization is working correctly, genes like actin and GAPDH — which are constitutively expressed at similar levels across all cell types and conditions — should show consistent medians. A systematic shift between Regen and NonRegen would indicate a normalization failure, not biology. ACTB is absent from the axolotl annotation; ACTG1 (gamma-actin, highly conserved) is the proxy.

---

### Slide 9: Section divider — Dataset 2: Xenopus Lin 2021
**Title**: Dataset 2: Xenopus laevis (Lin et al. 2021)

**Type**: section

---

### Slide 10: Raw 10x data goes through a four-step pipeline before any comparison
**Title**: Xenopus pipeline: four steps convert raw reads into a comparison-ready pseudobulk matrix

**Content**:
- **Step 1 — QC**: remove empty droplets (too few genes), doublets (too many genes), dying cells (high mito %) before any clustering — otherwise noise drives the biology
- **Step 2 — Normalization**: log-normalize so that expression reflects biology, not sequencing depth variation between cells
- **Step 3 — Harmony**: correct for library-level batch effects across samples without erasing biological differences between timepoints
- **Step 4 — Homeolog summing**: Xenopus laevis is allotetraploid — ptch1.L and ptch1.S are both expressed; sum them into one PTCH1 value before any cross-species comparison

**Notes**: Each step solves a specific technical problem. Skipping any one of them would introduce a systematic artifact: skip QC and empty droplets form a spurious "ghost" cluster; skip normalization and high-depth cells dominate; skip Harmony and libraries separate by sequencing run not cell type; skip homeolog summing and you undercount Xenopus expression relative to species with one gene copy.

---

### Slide 11: QC filtering removes technical artifacts before any biology is analyzed
**Title**: QC step: droplet capture is imperfect — empty droplets and doublets must be removed first

**Figure**: [Xenopus pre-filter QC violin](figures/xenopus_qc_violin.png)

**Type**: figure

**Figure Caption**: Xenopus (Lin et al. 2021) — pre-filter QC. nFeature > 500, nCount > 1,000, mito < 20%.

**Notes**: The violin shows the full distribution before filtering — this is why the thresholds are set where they are. Empty droplets contain ambient RNA, not a real cell — they appear as the peak at very low nFeature/nCount. Doublets (two cells captured together) appear as the tail at very high nFeature — they would cluster as artificial intermediate cell types if not removed. High mito% indicates a cell whose cytoplasm leaked out, leaving only mitochondria — these cells are stressed or dead and their transcriptomes are not informative.

---

### Slide 12: Harmony removes sequencing batch effects without erasing biology
**Title**: Harmony QC: timepoints interleave within clusters — batch corrected, biology intact

**Figure**: [Xenopus UMAP by timepoint](figures/xenopus_umap_timepoint.png)

**Figure Position**: right

**Figure Caption**: Post-Harmony UMAP colored by timepoint. Interleaving confirms batch correction worked.

**Content**:
- **Problem**: each 10x run introduces technical variation — without correction, cells cluster by run, not cell type
- **Harmony**: adjusts PCA embeddings so cells of the same type align across runs
- **QC**: every cluster must contain cells from multiple timepoints — single-timepoint clusters are batch artifacts

**Notes**: The intuitive check is: look at a cluster and see if it contains cells from 0dpa, 3dpa, 7-14dpa etc. If it does, Harmony worked — that cluster is a real cell type, not a sequencing artifact. If a cluster is monochromatic (one timepoint only), it is likely a batch-specific artifact. The timepoint UMAP is the QC artifact we save to verify this every run.

---

### Slide 13: BL and LBst cells separate in UMAP — confirming two biologically distinct programs
**Title**: Froglet (BL) and tadpole (LBst) occupy distinct UMAP regions despite being processed together

**Figure**: [Xenopus UMAP by condition](figures/xenopus_umap_condition.png)

**Figure Position**: right

**Content**:
- BL (froglet, non-regen) and LBst (tadpole, regen) separate in UMAP despite being processed together
- Some types shared (CT_progenitor); others exclusive (CT_fibroblast_E = BL only)
- Separation is biological — Harmony corrected within-condition batch without collapsing the between-condition biology

**Notes**: This is an important design choice — processing BL and LBst together with Harmony allows us to see whether their transcriptomes overlap or separate. If they completely mixed, the two conditions would be indistinguishable. Their separation here tells us the two biological programs are meaningfully different and justifies treating them as separate columns in the correlation matrix.

---

### Slide 14: Cell type annotation requires human expertise — canonical markers are the ground truth
**Title**: Cluster annotation uses canonical marker genes, not automated labels — here is the audit trail

**Figure**: [Xenopus canonical marker dotplot](figures/xenopus_marker_dotplot.png)

**Type**: figure

**Figure Caption**: Xenopus marker dotplot used to verify cluster-to-cell-type assignments (Lin et al. 2021).

**Notes**: This step is conceptually critical and often underappreciated. Unsupervised clustering groups cells by transcriptional similarity — it does not know what cell types those groups represent. Assigning labels requires looking at which genes are highest in each cluster and matching to known biology. The dotplot is saved every pipeline run because it is the only record of why a cluster was labeled the way it was. If a label seems wrong later, this dotplot is what you interrogate. Automated annotation tools (SingleR, Azimuth) exist but are less reliable for non-human, non-mouse species.

---

### Slide 15: Section divider — Dataset 3: Mouse digit tip
**Title**: Dataset 3: Mouse digit tip (Johnson et al. 2020)

**Type**: section

---

### Slide 16: Mouse pipeline mirrors Xenopus — same four steps, adjusted thresholds
**Title**: Mouse pipeline: same QC → normalization → Harmony → clusters framework as Xenopus

**Content**:
- **Input**: 10x CellRanger output, GSE143888 — 5 timepoints (0, 11, 12, 14, 17 dpa), regenerating digit tips
- **QC**: nFeature > 200, nCount > 500, mito < 25% — thresholds are lower because digit tip cells have fewer UMI counts than limb tissue; using Xenopus thresholds would discard real cells
- **Harmony**: aligns the 5 timepoint libraries; same logic as Xenopus — timepoints should interleave within clusters
- **Ortholog strategy**: toupper() only — "Ptch1" → "PTCH1"; works for ~95% of mouse-human gene pairs because most symbols differ only in capitalization

**Notes**: The QC thresholds differ from Xenopus for a biological reason — digit tip dissociations yield fewer transcripts per cell than whole limb preps. Using the same thresholds as Xenopus would remove a large fraction of real cells. This is why QC parameters are always dataset-specific, never copy-pasted.

---

### Slide 17: QC and Harmony produce a clean mouse dataset across all five timepoints
**Title**: Mouse digit tip: QC passes and 37,743 cells are retained across all five dpa timepoints

**Figure**: [Mouse pre-filter QC violin](figures/mouse_qc_violin.png)

**Figure Position**: right

**Figure Caption**: Mouse digit tip (GSE143888) — pre-filter QC. nFeature > 200, nCount > 500, mito < 25%.

**Content**:
- 37,743 cells retained after filtering — consistent with published cell counts
- Clusters not yet annotated: shown by number pending canonical marker review
- Same Harmony QC logic applies: timepoints must interleave within clusters

**Notes**: The mouse dataset is not yet annotated with cell type labels. Clusters are shown as numbers (C0–C18). The cell type annotation step is next in the pipeline — it requires inspecting canonical markers for mouse connective tissue, bone, immune, and vascular cell types.

---

### Slide 18: Harmony-integrated mouse UMAP: all five timepoints mix within every major cluster
**Title**: All five mouse timepoints contribute to every cluster — Harmony integration is successful

**Figure**: [Mouse UMAP by timepoint](figures/mouse_umap_timepoint.png)

**Figure Position**: right

**Content**:
- Cells from all 5 timepoints (0, 11, 12, 14, 17 dpa) are present within every major cluster
- No cluster is dominated by a single timepoint — same Harmony QC check as Xenopus
- Once annotation is complete, this becomes a cell type UMAP with timepoint overlays per cell type

**Notes**: The mouse UMAP is currently cluster-labeled because annotation is pending. The important QC message is that Harmony worked — the integration is technically sound, so when we do annotate we can trust that clusters represent biology, not sequencing batches.

---

### Slide 19: Section divider — Cross-species gene harmonization
**Title**: Cross-species gene name harmonization

**Type**: section

---

### Slide 20: Each species uses different gene names — a mapping strategy is required for comparison
**Title**: Without gene name harmonization, the three species share zero genes — three strategies solve this

**Content**:
- **Why this matters**: gene names are species-specific conventions; "Ptch1" (mouse), "ptch1.L" (Xenopus), and "PTCH1" (human/axolotl collab) refer to the same gene but would never intersect by string matching without transformation
- **Axolotl collab → HGNC**: already named — no mapping needed; 29,653 HGNC genes available
- **Xenopus → HGNC**: strip `.L`/`.S` suffix + toupper; both subgenome copies summed first. 36,487 homeolog rows → 30,134 unique base symbols
- **Mouse → HGNC**: toupper only — "Ptch1" → "PTCH1"; works for ~95% of pairs; genuinely divergent genes fall out at intersection

**Notes**: This is the step that makes the whole comparison possible. The limitation of the toupper approach for mouse is that some genes have genuinely different names between mouse and human — e.g., the mouse Ly6 family genes have no direct human counterpart by name. A biomaRt upgrade using the Ensembl ortholog database would recover these. For the GOI and HH genes in this project, all three strategies produce correct results because these well-conserved genes use consistent naming across species.

---

### Slide 21: A common gene set makes rho values directly comparable across all three species
**Title**: Using the same genes for all species removes a statistical confound from the correlation

**Content**:
- **The problem**: per-pairwise gene counts differ (Xenopus: 9,895 genes; Mouse: 10,987; Axolotl: 9,254) — a species with more shared genes could get a higher rho by chance, even if less biologically similar
- **Solution**: compute a second correlation using only genes shared across ALL datasets simultaneously (4,545 genes) — every species uses the same gene universe
- **Validation**: per-pairwise and common-set rankings agree — the Xenopus 7-14 dpa match is the top hit in both, confirming the result is real and not a gene-count artifact

**Notes**: This is a deliberate design choice, not just a nice-to-have. If the two strategies had disagreed — if Xenopus was top hit per-pairwise but axolotl was top hit in the common set — we would have to conclude that the per-pairwise result was an artifact of having more matching genes. Agreement between strategies is the validation. The common gene set is more conservative (fewer genes = less statistical power) but more interpretable.

---

### Slide 22: Section divider — Pseudobulk and correlation
**Title**: Pseudobulk construction and Spearman correlation

**Type**: section

---

### Slide 23: Pseudobulk averaging makes single-cell data comparable to bulk RNA-seq
**Title**: Pseudobulk: averaging all cells per timepoint converts sparse single-cell data into a bulk-like profile

**Content**:
- **Why**: individual cells are too sparse to correlate with bulk — most genes are 0 in any single cell due to dropout; you cannot compare a bulk profile to a single cell profile
- **How**: for each timepoint, compute the mean log-normalized expression across all cells of that timepoint → a dense gene × timepoint matrix that matches bulk RNA-seq format
- **Why mean across all cell types**: we are asking "which overall transcriptional state does BMS resemble?" — not which specific cell type drives it (that is a downstream question)

**Notes**: Pseudobulk is conceptually the inverse of what most people think of scRNA-seq for — we are collapsing the single-cell resolution to get a bulk-style profile. The point is that we still benefit from the scRNA-seq cell type annotation for downstream analyses; here we are using the dataset as a state classifier. An alternative would be to pseudobulk per cell type and then correlate, which we can do as a next step.

---

### Slide 24: Spearman rank correlation handles scale differences that would break Pearson
**Title**: Spearman rank correlation is used because expression scales differ fundamentally across species

**Content**:
- **Problem**: bulk RNA-seq normalized counts, scRNA-seq log-normalized values, and expression levels from different species are not on comparable scales — Pearson correlation would be dominated by whichever genes happen to be highest in each dataset
- **Spearman solution**: replace each gene's expression value with its rank within the profile — now every gene is on a 1-to-N scale regardless of species or assay, and the correlation measures whether the same genes tend to be relatively high (or low) in both profiles
- Result: a bulk conditions × scRNA-seq timepoints matrix of rho values, one per gene-set strategy

**Notes**: A concrete example of why Pearson fails here: ribosomal protein genes are always among the highest expressed in any dataset. If two profiles both have ribosomal genes at the top, Pearson gives a high correlation driven entirely by those genes — which tells you nothing biologically interesting. Spearman's ranks mean that if ribosomal genes are top-ranked in both profiles, they still contribute, but so does every other gene in proportion to its relative rank. The result is a more holistic measure of transcriptional state similarity.

---

### Slide 25: The heatmap shows consistent signal — all bulk conditions match Xenopus froglet 7-14 dpa
**Title**: All four bulk conditions most strongly resemble Xenopus froglet 7-14 dpa

**Type**: figure

**Figure**: [Bulk vs scRNA-seq correlation heatmap](figures/heatmap_perpairwise.png)

**Figure Caption**: Spearman rho — bulk conditions (rows) vs scRNA-seq pseudobulk (columns), per-pairwise genes.

**Notes**: Walk through the axes. Rows = 4 BMS conditions (top block) and 4 EtOH conditions (bottom block). Columns = all scRNA-seq timepoints grouped by species. The warm tile at LinBL_7-14dpa is consistent across all 8 bulk conditions — BMS and EtOH both match the same Xenopus state, which suggests the match reflects shared wound healing biology rather than the specific BMS treatment.

---

### Slide 26: BMS matches non-regenerative wound closure — not a regenerative blastema state
**Title**: BMS bulk data matches Xenopus non-regenerative wound closure (rho = 0.640), not blastema

**Figure**: [Reamputation best match per species](figures/reamputation_bestmatch.png)

**Figure Position**: right

**Figure Caption**: Best-matching timepoint per species for the Reamputation bulk condition.

**Content**:
- Top match: Xenopus froglet 7-14 dpa (rho = 0.640) — a non-regenerative wound-healing state
- Axolotl blastema and regenerative tadpole score lower (rho 0.47-0.56)
- BMS inhibits NF-kB at a window when the axolotl transcriptome resembles wound healing, not blastema

**Notes**: Frame carefully — this is an initial finding from the pipeline validation, not a final conclusion. The observation is reproducible and biologically coherent (Xenopus froglet 7-14 dpa has been characterized as a fibrotic wound healing response). What it means for BMS biology requires more analysis — the GOI and HH pathway comparisons are the next step.

---

### Slide 27: Every step of the pipeline has a QC checkpoint and produces an inspectable artifact
**Title**: The pipeline is reproducible, validated at each step, and extensible to new datasets

**Content**:
- **Audit trail**: every step saves an inspectable figure (QC violin, UMAP by sample, marker dotplot)
- **Logged gene counts**: gene intersections and rho values logged per run — no silent failures
- **Extensible**: new datasets plug in via one `pseudobulk_*()` function; rest of pipeline unchanged
- **Next**: mouse annotation, biomaRt ortholog upgrade, GSE135985 non-regen mouse, GOI/HH analyses

**Notes**: The pipeline is the lab's contribution at this stage — not just this one result. The emphasis on QC checkpoints means that when results change (new data, new annotation), the reason for the change can be traced back to a specific step.

---
