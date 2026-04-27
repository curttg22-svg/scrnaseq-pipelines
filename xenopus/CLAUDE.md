# Xenopus Limb Regeneration scRNA-seq Pipeline

## Overview
Five-timepoint Xenopus laevis limb regeneration dataset (0, 3, 7-14, 14, 14-52 dpa).
Two GEO batches (GSM5045xxx + GSM5057xxx) integrated with Harmony by sample.

## Annotation Strategy
Cell type labels are assigned by **marker-based verification**, not by assumed cluster order.

**Workflow:**
1. Run `01_xenopus_qc_clustering.R` to produce clustered object with Harmony embedding
2. Run `02_xenopus_annotation.R` — this generates `xenopus_canonical_marker_dotplot.pdf` first
3. **Inspect the dotplot before accepting labels.** Each cluster's dominant markers must match the assigned cell type. If they do not, correct the `RenameIdents` call and re-run.
4. Re-run GOI visualization (`06_goi_visualization.R`) after any annotation change

**Key marker panel** (canonical_markers in 02_xenopus_annotation.R):
- col1a1/vim → CT fibroblast
- col2a1/sox9 → cartilage/chondrocyte
- krt8/krt18/epcam → epithelial/keratinocyte
- hba1/hbb → erythrocyte (included to resolve ambiguous clusters)
- ptprc → pan-immune
- pecam1/cdh5 → endothelial
- acta2/tagln → smooth muscle
- sox10 → neural crest / Schwann

## Clusters to Re-verify After Re-clustering
After any re-clustering (resolution change, new samples, new Harmony run):
- Cluster labeled Erythrocyte: confirm hba1/hbb expression; col1a1 signal in this row indicates possible mislabel
- Cluster labeled CT_cartilage_assoc: confirm col2a1/sox9; epcam signal indicates possible mislabel
- Any cluster with unexpected sox10: reassign to Neural/Schwann category

## Workflow Conventions
- Do not rewrite existing scripts; extend or modify in place
- Verify prior outputs before starting new analysis
- Use ASCII hyphens in plot titles (not em-dashes) to avoid R encoding warnings
- All commits go to personal repo only (curttg22-svg/scrnaseq-pipelines); course is complete
