# Cross-Species Regeneration Comparison Pipeline

## Workflow Guidelines

When extending existing analyses (PCA, atlas annotation, etc.), read and reuse the existing script rather than writing from scratch. Verify prior feedback/TODOs have been addressed before adding new work.

## Bioinformatics Workflow

Before starting large downloads (>1 GB) or installing new bioinformatics tools (Salmon, BWA, etc.), confirm the plan with the user first. Check if data/tools already exist locally before downloading.

### Pre-flight Protocol

Before any operation that downloads >1 GB, installs new tools, or runs >10 minutes, output a structured plan and wait for approval:

**(a) Inputs and sources** — what files are needed and where they come from  
**(b) Disk space delta** — how much space will be used at peak and what remains after cleanup  
**(c) Wall-clock estimate** — realistic time to complete  
**(d) Existing assets check** — search the repo first; confirm no local file, script, or cached result already covers this  
**(e) Minimal vs. comprehensive approach** — state the lightest path that achieves the goal, and what the fuller approach would add

Do not execute until the user approves.

## R / Seurat Conventions

For Seurat v5 objects, use `unname()` when accessing assay data, and always verify that columns referenced in re-plot scripts are actually saved in the RDS file before running.

## Plotting Conventions

Avoid em-dashes (—) in plot titles, labels, and other rendered text; use regular hyphens or colons instead to prevent encoding warnings.
