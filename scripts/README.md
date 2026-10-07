# Cleaned and annotated scripts

These are the main working scripts. `archive/` preserves the supplied originals for provenance. Each script has explicit inputs, comments explaining the analysis, and a source/adaptation record below. The numbered documents in `analyses/` explain how they connect to the paper.

The scripts are prepared for execution with the required data and packages. Scientific results have not been rerun or matched to final paper outputs. Historical settings below reconstruct a supplied branch; they are not asserted to be the final publication settings.

## Script index

| Script | Purpose | Source and treatment |
| --- | --- | --- |
| [alignment/starsolo.sh](alignment/starsolo.sh) | STARsolo counts | Original starsolo.sh settings retained; paths, threads and zero-based sample index are arguments; quoting, validation and output checks added. Removed inactive experiments and personal scheduler settings. |
| [snrna/01_import_qc_sct.R](snrna/01_import_qc_sct.R) | Matrix import, metadata, merge, QC and SCT | Reconstructs the COMMENTED March upstream branch as explicit stages. QC and regression settings are required arguments. Sample metadata is supplied rather than inferred. |
| [snrna/02_integrate.R](snrna/02_integrate.R) | SCT integration, PCA, UMAP, clustering and save | March active integration tail; optional frozen barcode selection also makes the recovered SOP lineage branch usable. Fixed missing parenthesis in lineage FindNeighbors by using a complete shared integration step. |
| [snrna/03_cell_cycle.R](snrna/03_cell_cycle.R) | S/G2M scoring | Clean scoring block. Exports marker audit, explicitly excludes unmatched IDs and deduplicates within phase; original mapping remains unchanged. This effective-feature treatment needs comparison with the historical run. |
| [snrna/04_classify_ivory.R](snrna/04_classify_ivory.R) | RNA ivory-high/low assignment | Complete race-specific threshold branch; assignment occurs after defining the selected object. Thresholds must be supplied. Does not invent the incomplete RNA_by_Ivory threshold. |
| [snrna/05_find_markers.R](snrna/05_find_markers.R) | DE table | Clean nucleus-level FindMarkers branch with explicit assay, both groups, filter settings and cell membership. Explicit Wilcoxon replaces the omitted default; does not claim to implement pseudobulk. Resolves the Markers/Markers2 name collision by writing one named output. |
| [snrna/06_expression_ranking.R](snrna/06_expression_ranking.R) | S8 aggregate CPM ranking | Clean AggregateExpression/edgeR block; selected input controls cohort. Exports ranked values and a vector PDF. |
| [plots/lineage_panels.R](plots/lineage_panels.R) | UMAP, feature and proportion panels | Rebuilds the plot types using explicit metadata fields and frozen cell selection. Proportions are among all selected nuclei per group; a category filter must be applied explicitly upstream. Exports membership and coordinates. |
| [plots/volcano.R](plots/volcano.R) | DE volcano PDF | Rebuilds plot from saved DE table; explicit thresholds, optional labels/categories, adjusted p-value floor only for displaying numerical zeros. |
| [plots/dotplot.R](plots/dotplot.R) | DE gene dot plot | Explicit gene/label table and group membership; Seurat feature IDs are not renamed. Exports plotted values. |
| [plots/go_enrichment.R](plots/go_enrichment.R) | DirectGO result PDF | Clean CSV plot with explicit correction columns, threshold and title; enrichment statistics must come from the original service exports. |
| [locus/extract_spliced_reads.sh](locus/extract_spliced_reads.sh) | Strand-specific spliced BAMs | Cleans the active samtools pipeline only; input/output checks, threading and temporary files added. Earlier commented locus extraction and barcode selection are not guessed. |

## Input formats

Use full Seurat cell names, including sample prefixes, consistently across membership files. RDS and `.h5seurat` input objects are supported. R scripts save options and `sessionInfo()` next to outputs. Most scripts refuse to overwrite the principal output; use a new output path for each parameter branch. Stochastic steps require an explicit seed, newly added for reproducibility because none was recorded in the source.

`sample_id`, `counts_dir`, `time`, `race`, `genotype` are the required columns in a sample TSV. Paths are interpreted relative to the current working directory. A header-only template is in [config/samples.template.tsv](../config/samples.template.tsv). Populate it from the confirmed sample sheet. Do not silently substitute older developmental times.

Barcode selection files have one full cell name per line, with no header. Contrast membership TSVs have `barcode` and `group` columns. Gene-label TSVs have `gene_id` and `label` columns. Optional volcano annotation TSVs have `gene_id` and `category` columns; they must be derived from the confirmed target/orthology/locus tables.

## Example execution: historical March reconstruction

Run from the repository root in the appropriate environment. Replace the example input paths with actual files. These commands select the commented 6,000-variable-feature/time-regressed upstream branch and the active 6,000-feature/four-PC integration tail. Confirm that pairing against the final object/logs before declaring a reproduced paper analysis.

```bash
Rscript scripts/snrna/01_import_qc_sct.R \
  --samples inputs/samples.tsv \
  --ribo inputs/Hmel_ribosomal_genes.txt \
  --mito inputs/Hmel_mito_genenames.txt \
  --output results/Hm_merged.rds \
  --min-cells 3 --min-features 2000 --qc-features 200 --max-mt 20 \
  --variable-features 6000 --regress time --seed 1

Rscript scripts/snrna/02_integrate.R \
  --input results/Hm_merged.rds --output results/Hm_combined.sct.rds \
  --split-by race --features 6000 --pcs 50 --dims 4 --resolution 1 --seed 1 \
  --h5seurat results/Hm_combined.sct.h5seurat
```

To use an already prepared object, start at integration. A frozen SOP barcode selection can be supplied with `--cells inputs/sop_cells.txt`. The recovered lineage trial used `--features 7000 --dims 3 --resolution 1`; this is a separate historical branch, not confirmation of final settings.

```bash
Rscript scripts/snrna/03_cell_cycle.R \
  --input results/Hm_combined.sct.rds \
  --genes archive/original/ivory_proj_transfer_aug2025/CellCycleGenes.txt \
  --assay RNA --output results/Hm_cell_cycle.rds

Rscript scripts/snrna/04_classify_ivory.R \
  --input inputs/annotated_lineage.rds --cells inputs/early_scale_cells.txt \
  --assay RNA --iv-threshold 1 --other-threshold 4 \
  --output results/early_classified.rds

Rscript scripts/snrna/05_find_markers.R \
  --input results/early_classified.rds --membership inputs/early_hi_vs_lo.tsv \
  --assay RNA --group1 hi --group2 lo --min-pct 0.01 --logfc-threshold 0.1 \
  --output results/early_hi_vs_lo.tsv
```

The DE filter values above are explicit EXAMPLE settings; recover the executed Seurat version/filter values before matching results. For WT versus mutant, supply membership limited to confirmed ivory-high nuclei and use `--group1 wt --group2 mut`. Run late contrasts with their separate memberships. The script does not infer genotype or contrast inclusion.

```bash
Rscript scripts/plots/volcano.R \
  --input results/early_hi_vs_lo.tsv --p-cutoff 0.05 --fc-cutoff 1 \
  --output results/early_hi_vs_lo.volcano.pdf

Rscript scripts/plots/dotplot.R \
  --input results/early_classified.rds --membership inputs/early_hi_vs_lo.tsv \
  --genes inputs/early_gene_labels.tsv --assay RNA \
  --output results/early_hi_vs_lo.dotplot.pdf

Rscript scripts/plots/go_enrichment.R \
  --input inputs/UpGenes_DirectGoEnrichment_PANGEA.csv \
  --filter-column Benjamini...Yekutieli --color-column Benjamini...Yekutieli \
  --cutoff 0.05 --title 'Early upregulated genes' --output results/early_up.go.pdf

Rscript scripts/snrna/06_expression_ranking.R \
  --input inputs/final_selected_lineage.rds --group-by celltypes --group scale \
  --output results/scale_expression_ranking.tsv

Rscript scripts/plots/lineage_panels.R \
  --input inputs/annotated_lineage.rds --assay RNA --group-by race \
  --category celltypes --feature Ivory --seed 1 --output-dir results/lineage_panels
```

Figure palettes, dimensions and annotation categories in rebuilt plot scripts are initial clean equivalents, not exact final artwork. Verify them against the final vector PDFs. High/low thresholds, final cohort, enrichment correction and DE statistical description remain explicit unresolved scientific choices.

## Alignment and splice extraction

```bash
bash scripts/alignment/starsolo.sh \
  inputs/jobs 0 inputs/fastqs inputs/STAR_index inputs/edited.gtf results/STAR 64

bash scripts/locus/extract_spliced_reads.sh \
  inputs/iv_20.CortexLocus.fw.bam inputs/iv_20.CortexLocus.rv.bam \
  results/iv_20.CortexLocusSplice 12
```

STAR prefixes now end in a dot (`redo.<sample>.`) to separate names cleanly; document this naming change when locating matrices. Counting parameters and R2/R1 order are unchanged. The original inactive BAM/barcode options remain inactive; the complete per-barcode coverage workflow needs the actual BAM-producing command and missing coverage scripts.

## Software and validation

Requirements: Bash, awk, STAR, samtools; R with Seurat, sctransform, SeuratDisk for H5 input/export, edgeR for CPM, ggplot2 and ggrepel for labelled volcano plots. Exact versions still need the original environment. Multi-layer Seurat inputs may need a verified layer preparation step; scripts do not silently join layers. The historical 2023 glmGamPoi/cell-cycle-regression experiments remain archived and have not been collapsed into this March branch.

See [CLEANUP_NOTES.md](../docs/CLEANUP_NOTES.md) for changes, checks and current limitations. No replacement RNAhybrid, Nanopore or reference-editing algorithms have been invented for missing code.
# Drive additions (7 October 2026)

`locus/plot_barcode_reads.R` cleans the recovered read/ribbon loop. It requires Rsamtools, GenomicRanges, IRanges, ggplot2 and scico. Run from the repository root with an indexed strand-filtered BAM containing CB tags:

```bash
Rscript --vanilla scripts/locus/plot_barcode_reads.R \
  --bam /path/to/iv_20.CortexLocus.fw.bam \
  --metadata data/drive/iv_reanalysis/bcdf.all.txt \
  --sample iv_20 --contig Hmel215003o \
  --start 1400000 --end 1600000 --output-prefix results/locus/iv_20
```

Set the interval to match the actual extraction. Outputs are vector read/ribbon PDFs, source tables and matching counts. Ribbons represent the historical occupied-position calculation, not basewise coverage; see [DRIVE_INTAKE.md](../docs/DRIVE_INTAKE.md).

To recover tables and metadata from a saved workspace:

```bash
Rscript --vanilla tools/inspect_workspace.R /path/to/workspace.RData workspace_inventory
```

This creates an object inventory and table/Seurat metadata exports without sourcing analysis scripts. Saved Seurat objects require compatible packages and enough memory. This utility has not yet been executed against a project workspace.


## Compact saved-workspace recovery

The `.RData` is now available in Drive but is too large for this runtime. In an R session with the workspace already loaded, source `tools/export_workspace_audit.R` and call `export_workspace_audit(outdir = "ivory_workspace_audit")`. This avoids another full workspace load and exports tables/metadata without matrices. The fresh-process `tools/inspect_workspace.R` wrapper remains available on a sufficiently provisioned machine. See [WORKSPACE_AUDIT.md](../docs/WORKSPACE_AUDIT.md) for findings and exact export instructions.

For the recovered late GO PDF candidate, history and panel term labels support the `_3` files with BY filtering for both panels, BY coloring for up and BH coloring for down:

```bash
Rscript scripts/plots/go_enrichment.R \
  --input data/drive/iv_reanalysis/UpGenesLate_DirectGoEnrichment_PANGEA_3.csv \
  --filter-column Benjamini...Yekutieli --color-column Benjamini...Yekutieli \
  --cutoff 0.05 --title 'Late upregulated genes' --output results/late_up.go.pdf

Rscript scripts/plots/go_enrichment.R \
  --input data/drive/iv_reanalysis/DownGenesLate_DirectGoEnrichment_PANGEA_3.csv \
  --filter-column Benjamini...Yekutieli --color-column Benjamini...Hochberg \
  --cutoff 0.05 --title 'Late downregulated genes' --output results/late_down.go.pdf
```

These produce separate clean PDFs, not the original cowplot assembly. Preserve the 27 up rows even though there are only 23 distinct display labels. This matching supports a saved panel candidate, not a complete final-paper validation.
