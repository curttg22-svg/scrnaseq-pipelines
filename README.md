# scrnaseq-pipelines

Cross-species single-cell RNA-seq analysis of limb regeneration in *Xenopus laevis*, *Ambystoma mexicanum* (axolotl), and mouse digit tip. Developed at Wake Forest University.

The pipeline proceeds in two stages: (1) per-species QC, clustering, and annotation, then (2) cross-species integration, ortholog mapping, pseudobulk construction, and dimensionality reduction.

---

## Repository Structure

```
scrnaseq-pipelines/
├── xenopus/              # X. laevis — Lin et al. 2021 (GSE165901)
├── mouse-digit/          # Mouse digit tip — Johnson et al. 2020 (GSE143888)
├── axolotl-leigh-2018/   # A. mexicanum — Leigh et al. 2018 (SCP422/489/499/500)
├── axolotl-gerber-2021/  # A. mexicanum — Gerber et al. 2021 (GSE165901)
├── comparison/           # Cross-species atlas, ortholog mapping, pseudobulk MDS
└── reference/            # Shared ortholog tables and annotation helpers
```

---

## Part I -- Individual Species Pipelines

Run each pipeline from within its folder. Scripts are numbered in execution order.

---

### 1. Xenopus laevis -- Lin et al. 2021

**Source:** GEO [GSE165901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165901) -- *Dev. Cell* 2021  
**Biology:** Froglet (non-regenerative) and tadpole NF50-52 (regenerative) limb tissue

#### Download

```r
# xenopus/
source("R/00_download_data.R")
```

Downloads all samples from GEO using `GEOquery`. Two experimental groups:

| Sample | Stage | Regenerative? |
|--------|-------|---------------|
| BL 0 dpa, 3 dpa, 7-14 dpa (pool), 14 dpa, 14-52 dpa (pool) | Froglet | No |
| LBst NF50, NF51, NF52 | Tadpole | Yes |

Excluded: GSM5057666 (mixed regen + non-regen pool), GSM5057667 (transplant experiment).

#### QC, Normalization, Clustering, UMAP

```r
source("R/01_xenopus_qc_clustering.R")
```

- Per-sample filtering on `nFeature_RNA`, `nCount_RNA`, and `percent.mt`
- `NormalizeData` -> `FindVariableFeatures` (3,000 HVGs) -> `ScaleData` -> `RunPCA`
- Harmony integration across samples (`RunHarmony(group.by.vars = "sample")`)
- `RunUMAP` on Harmony embeddings -> `FindNeighbors` -> `FindClusters`
- Output: `results/xen_merged.rds`

#### Cell Type Annotation and Dot Plots

```r
source("R/02_xenopus_annotation.R")
```

Generates a canonical marker dot plot first -- **inspect before accepting labels.** Each cluster's dominant markers must match the assigned cell type.

| Marker | Cell type |
|--------|-----------|
| `col1a1`, `vim` | CT fibroblast |
| `col2a1`, `sox9` | Cartilage / chondrocyte |
| `krt8`, `krt18`, `epcam` | Epithelial / keratinocyte |
| `hba1`, `hbb` | Erythrocyte |
| `ptprc` | Pan-immune |
| `pecam1`, `cdh5` | Endothelial |
| `acta2`, `tagln` | Smooth muscle |
| `sox10` | Neural crest / Schwann |

Output: `results/xen_merged.rds` (with `cell_type` labels), `results/xenopus_canonical_marker_dotplot.pdf`

#### Additional Xenopus Scripts

| Script | What it does |
|--------|-------------|
| `02b_xenopus_early_timepoints.R` | Sub-clusters early tadpole timepoints |
| `06_goi_visualization.R` | UMAP feature plots and violin plots for GOIs |
| `07_homeolog_expression.R` | L vs. S subgenome homeolog expression comparison |

---

### 2. Mouse Digit Tip -- Johnson, Masias & Lehoczky 2020

**Source:** GEO [GSE143888](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143888) -- *Dev. Cell* 2020  
**Biology:** Adult mouse digit tip regeneration (distal phalanx amputation), 5 timepoints

#### Download

```r
# mouse-digit/
source("R/00_download_data.R")
```

Downloads 5 10x Genomics CellRanger output directories from GEO:

| GSM | Timepoint |
|-----|-----------|
| GSM4276219 | 0 dpa (unamputated control) |
| GSM4276220 | 11 dpa |
| GSM4276221 | 12 dpa |
| GSM4276222 | 14 dpa |
| GSM4276223 | 17 dpa |

#### QC, Normalization, Clustering, UMAP

```r
source("R/01_mouse_qc_clustering.R")
```

- Per-sample QC filtering, Harmony integration by `sample`
- Output: `results/mouse_merged.rds`

---

### 3. Axolotl -- Leigh et al. 2018

**Source:** Broad Single Cell Portal -- [SCP422](https://singlecell.broadinstitute.org/single_cell/study/SCP422) / [SCP489](https://singlecell.broadinstitute.org/single_cell/study/SCP489) / [SCP499](https://singlecell.broadinstitute.org/single_cell/study/SCP499) / [SCP500](https://singlecell.broadinstitute.org/single_cell/study/SCP500)  
Also mirrored at GEO [GSE121737](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121737)  
**Publication:** Leigh et al. 2018, *Nat. Commun.* doi:10.1038/s41467-018-07604-0  
**Biology:** Axolotl forelimb across four stages -- Intact, Wound Healing, Early Blastema, Mid Blastema

#### Download

Download each SCP study from the Broad portal. Place matrices in:

```
axolotl-leigh-2018/data/SCP422_intact/
axolotl-leigh-2018/data/SCP489_WH/
axolotl-leigh-2018/data/SCP499_EB/
axolotl-leigh-2018/data/SCP500_MB/
```

#### QC, Integration, Clustering, UMAP, Annotation

```r
# axolotl-leigh-2018/
source("R/01_load_scp_datasets.R")   # load SCP tab-separated format -> Seurat objects
source("R/02_integrate_harmony.R")   # Harmony integration across timepoints
source("R/03_annotation_umap.R")     # clustering, cell type labels, UMAP
```

Gene symbols in Leigh 2018 are Trinity transcriptome IDs mapped to SwissProt. The helper table `helper/trinity_to_symbol_map.csv` converts them to HGNC symbols automatically during cross-species steps.

---

## Part II -- Cross-Species Integration and PCA

All scripts below run from the `comparison/` directory.

---

### 4. Cross-Species Atlas (Harmony)

```r
# comparison/
source("R/build_cross_species_atlas.R")
```

Merges all three species into a shared atlas:

- **Axolotl collaborator dataset** (`data/axo_collab/`) -- 43K cells, regen + non-regen conditions, HGNC symbols. *Private/unpublished -- contact PI for access.*
- **Xenopus:** `xenopus/results/xen_merged.rds`
- **Mouse:** `mouse-digit/results/mouse_merged.rds`

Steps:
1. Subsample to 20,000 cells per species (configurable via `SUBSAMPLE_N`)
2. Normalize and find 3,000 HVGs per species
3. Merge and run Harmony (`group.by.vars = c("sample", "species")`)
4. Cluster at resolution 0.5 -> assign shared `cell_type_atlas` labels
5. Output: `results/cross_species_atlas.rds`

---

### 5. Ortholog Mapping

Gene names are harmonized to HGNC uppercase symbols via `R/utils_harmonize.R`:

| Species | Input format | Conversion |
|---------|-------------|------------|
| Axolotl (collab) | HGNC symbols | Pass-through |
| Xenopus | `col1a1.L` / `col1a1.S` homeologs | Strip `.L`/`.S`, `toupper()`, average homeologs |
| Mouse | Title-case (`Col1a1`) | `toupper()` |

The ~9,300 genes shared across all three species after harmonization form the **common gene set** used for pseudobulk and PCA.

---

### 6. Pseudobulk Construction

```r
source("bulk_vs_scrna_correlation.R")
```

Collapses each species' scRNA-seq to a **genes x timepoints pseudobulk matrix** (mean log-normalized expression per timepoint) and saves `results/pseudobulk_cache.rds`. This cache is the sole input to the MDS analysis.

---

### 7. Cross-Species Regenerative MDS

```r
source("R/plot_pca_exploration_noxentad.R")   # full exploration: all PC pairs x all normalizations
source("R/plot_regeneration_pca.R")            # main regenerative spectrum figure
```

Approach:
1. **Within-species centering** -- subtract per-species gene mean, removing species-identity signal so only regeneration-driven variation remains
2. **Spearman distance MDS** (cmdscale, k=6) on 12 reference pseudobulks (Xenopus tadpole excluded -- see rationale below)
3. **Projection** -- new samples placed into this space as Spearman-weighted centroids

**Why tadpoles are excluded from the reference:** Including NF50-52 tadpole samples inflates PC1 to 39.6% with species-identity variance, compressing the regenerative signal. With tadpoles excluded, variance redistributes evenly (PC1=24.2%, PC2=20.1%, PC3=18.0%) and each axis captures interpretable biology.

| Axis | Variance | Biology |
|------|----------|---------|
| PC1 | 24.2% | Regenerative gradient -- separates regen from non-regen (Wilcoxon p = 0.009) |
| PC2 | 20.1% | Metabolic demand vs. late maturation / ciliogenesis |
| PC3 | 18.0% | Proliferative blastema (negative) vs. fibrotic scar (positive) |

Output: `figures/pca_noxentad/`

---

## Dependencies

```r
# CRAN
install.packages(c("Seurat", "harmony", "ggplot2", "ggrepel", "dplyr",
                   "tidyr", "patchwork", "Matrix", "data.table", "stringr"))

# Bioconductor
BiocManager::install(c("GEOquery", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db"))
```

R >= 4.2, Seurat >= 5.0. Large files (`.rds`, raw count matrices) are excluded via `.gitignore` -- all results are reproducible from scripts alone.

---

## Notes

- Run all scripts from within their project folder so relative paths resolve correctly
- Seurat v5: use `unname()` when accessing assay data slots directly
- Xenopus gene names follow `.L`/`.S` homeolog suffix convention -- use `utils_harmonize.R`, do not strip manually
- Use hyphens in plot titles, not em-dashes, to avoid encoding warnings
- Before any download >1 GB or new tool install, follow the pre-flight protocol in `comparison/CLAUDE.md`
