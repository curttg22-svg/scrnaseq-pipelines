# scrnaseq-pipelines

Single-cell RNA-seq analysis pipelines for *Xenopus laevis* and *Ambystoma mexicanum*
(axolotl) limb regeneration. Developed at Wake Forest University.

## Datasets

| Folder | Species | Dataset | Accession |
|--------|---------|---------|-----------|
| `xenopus/` | *X. laevis* | Gerber et al. 2021 | GEO GSE165901 |
| `axolotl-gerber-2021/` | *A. mexicanum* | Gerber et al. 2021 (CT fate-map, 11 dpa) | GEO GSE165901 |
| `axolotl-leigh-2018/` | *A. mexicanum* | Leigh et al. 2018 (Intact → WH → EB → MB) | GEO GSE121737 / SCP422 |

## Workflow

Each dataset folder contains numbered R scripts meant to be run in order:

```
01_  Load raw data / QC
02_  Normalization, clustering, UMAP
03_  Cell type annotation
04+  Downstream analysis (GOI visualization, homeolog expression, LOESS curves)
```

Run scripts from within their dataset folder so relative paths resolve correctly.

## Dependencies

- R ≥ 4.2
- Seurat ≥ 5.0
- harmony, ggplot2, dplyr, patchwork, stringr, Matrix, data.table

## Notes

- Large files (`.rds`, raw data) are excluded from git via `.gitignore` — results
  are fully reproducible from the scripts
- Xenopus gene names follow the `.L` / `.S` homeolog suffix convention
- Axolotl (Leigh 2018) uses Trinity transcriptome IDs; gene symbols are extracted
  automatically via `axolotl-leigh-2018/helper/trinity_to_symbol_map.csv`
