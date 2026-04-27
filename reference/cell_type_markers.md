# Cross-Species Cell Type Marker Reference

Canonical marker genes used for cell type annotation across all scRNA-seq pipelines
in this repo. Human gene symbols are the shared key; species-specific gene IDs are
listed below each cell type.

---

## Gene Naming Conventions by Species

| Species | Gene ID format | Example | Source |
|---|---|---|---|
| Human (reference) | HGNC symbol | COL1A1 | Standard |
| *Xenopus laevis* | Symbol.L / Symbol.S (homeolog) | col1a1.L, col1a1.S | XenBase Xl10.1 |
| *A. mexicanum* (Gerber/cross-species) | AMEX60DD ID, mapped via bridge | AMEX60DD000123 | AmexT_v47-AmexG_v6.0-DD.gtf |
| *A. mexicanum* (Leigh 2018) | Trinity ID, mapped via helper CSV | c12345_g1_i1 | trinity_to_symbol_map.csv |

---

## Annotation Files

### Xenopus laevis
- **Genome assembly:** XenBase Xl10.1 (GCF_017654675.1)
- **Gene naming:** All genes carry .L (long/Chr1L) or .S (short/Chr1S) homeolog suffixes
- **Reference:** https://www.xenbase.org/xenbase/static-xenbase/ftpDatafiles.jsp
- **Pipeline script:** `xenopus/R/02_xenopus_annotation.R`

### Axolotl — cross-species pipeline (Gerber 2021 / cross-species comparison)
- **Transcriptome:** AmexT_v47, genome AmexG_v6.0-DD
- **GTF file:** `AmexT_v47-AmexG_v6.0-DD.gtf` (place at `data/annotation/` — not tracked in git, too large)
- **Download:** https://axolotlomics.org
- **Bridge built by:** `xenopus/R/04_axolotl_gtf_annotation.R` → `results/final_bridge.rds`
- **Format:** GTF gene_name field encodes `AXOLOTLNAME [NR]|HUMANNAME [HS]`; script extracts `[HS]` symbol
- **Pipeline script:** `xenopus/R/05_axolotl_reannotation.R`

### Axolotl — Leigh 2018 (multi-timepoint regeneration)
- **Transcriptome:** Custom Trinity assembly with embedded UniProt annotations
- **Format:** Trinity IDs contain `^sp|ACCESSION|SYMBOL_SPECIES^` patterns; Seurat sanitizes `|` to `-`
- **Bridge:** `axolotl-leigh-2018/helper/trinity_to_symbol_map.csv`
- **Symbol extraction strategy:** Prefer `_HUMAN` suffix; fall back to first species hit; GOI-specific aliases in `UNIPROT_ALIASES` dict
- **Pipeline script:** `axolotl-leigh-2018/R/01_load_scp_datasets.R`

---

## Canonical Marker Genes by Cell Type

These markers are used in canonical dotplots for annotation verification.
Always inspect the dotplot before accepting cluster labels — do not assume labels from cluster order.

### Connective Tissue / Fibroblast
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| COL1A1 | col1a1.L, col1a1.S | via final_bridge |
| COL1A2 | col1a2.L, col1a2.S | via final_bridge |
| VIM | vim.L | via final_bridge |
| FN1 | fn1.L | via final_bridge |
| DCN | dcn.L | via final_bridge |

### Cartilage / Chondrocyte
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| COL2A1 | col2a1.L | via final_bridge |
| SOX9 | sox9.L | via final_bridge |
| ACAN | acan.L | via final_bridge |

### Smooth Muscle / Myofibroblast
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| ACTA2 | acta2.L | via final_bridge |
| TAGLN | tagln.L | via final_bridge |
| MYH11 | myh11.L | via final_bridge |

### Epithelial / Keratinocyte
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| EPCAM | epcam.L | via final_bridge |
| KRT8 | krt8.L | via final_bridge |
| KRT18 | krt18.L | via final_bridge |
| KRT17 | krt17.L | via final_bridge |
| CDH1 | cdh1.L | via final_bridge |

### Macrophage
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| AIF1 (IBA1) | aif1.L | via final_bridge |
| CD68 | cd68.L | via final_bridge |
| PTPRC (CD45) | ptprc.L | via final_bridge |
| CSF1R | csf1r.L | via final_bridge |
| C1QA | c1qa.L | via final_bridge |

### Neutrophil
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| MPO (MPX) | mpx.L | via final_bridge |
| S100A8 | s100a8.L | via final_bridge |
| S100A9 | s100a9.L | via final_bridge |
| ELANE | elane.L | via final_bridge |

### B Cell
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| CD79A | cd79a.L | via final_bridge |
| MS4A1 (CD20) | ms4a1.L | via final_bridge |
| PAX5 | pax5.L | via final_bridge |

### T Cell
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| CD3E | cd3e.L | via final_bridge |
| CD3D | cd3d.L | via final_bridge |
| CD247 | cd247.L | via final_bridge |

### Endothelial
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| PECAM1 (CD31) | pecam1.L | via final_bridge |
| CDH5 | cdh5.L | via final_bridge |
| VWF | vwf.L | via final_bridge |
| CLDN5 | cldn5.L | via final_bridge |

### Erythrocyte
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| HBA1 | hba1.L | via final_bridge |
| HBB | hbb.S | via final_bridge |
| ALAS2 | alas2.L | via final_bridge |

> **Note:** hba1.L and hbb.S are included in the Xenopus canonical marker panel specifically
> to resolve ambiguity in clusters 14 and 16, which show unexpected col1a1/epcam signals.
> If either cluster lacks hba1/hbb expression, re-assign to chondrocyte or epithelial respectively.

### Neural Crest / Schwann Cell
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| SOX10 | sox10.L | via final_bridge |
| S100B | s100b.L | via final_bridge |
| NGFR | ngfr.L | via final_bridge |
| MPZ | mpz.L | via final_bridge |

### Proliferating (marker-positive cycling cells)
| Human | Xenopus (.L/.S) | Axolotl (via bridge) |
|---|---|---|
| MKI67 | mki67.L | via final_bridge |
| TOP2A | top2a.L | via final_bridge |
| PCNA | pcna.L | via final_bridge |

---

## Cell Type Harmonization Across Species

For cross-species GOI comparison, cell types are harmonized to shared categories:

| Harmonized label | Xenopus label(s) | Axolotl (Leigh 2018) label(s) | Axolotl (Gerber) label(s) |
|---|---|---|---|
| CT_fibroblast | CT_fibroblast_A, CT_fibroblast_B | SSC/fibroblast, fibroblast blastema | Blastema_CT_A/B, CT_differentiating |
| Smooth_muscle | Smooth_muscle | - | - |
| Chondrocyte | CT_cartilage_assoc | - | - |
| Epithelial | Keratinocyte, Keratinocyte_2/3, Basal_epithelial | Basal epidermis, intermediate epidermis | - |
| Macrophage | Macrophage, Macrophage_2 | Myeloid/macrophage, recruited macrophage | Macrophage |
| Neutrophil | Neutrophil | Neutrophil subtypes | - |
| Endothelial | Endothelial | Endothelial | - |
| Erythrocyte | Erythrocyte, Erythrocyte_2 | Erythrocyte/RBC | - |
| Neural/Schwann | Neural_crest, Neural_2, Schwann_cell, Neuron | Schwann | - |
| Proliferating | Proliferating_A, Proliferating_B | Proliferating epidermis | - |
| Immune_other | - | T cell, early B cell | - |
