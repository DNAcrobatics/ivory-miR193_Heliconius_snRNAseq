# Analysis scripts

Run commands from the repository root. Each script accepts explicit input paths and analysis settings. Preserve separate output directories for distinct cohorts, assays or parameter choices.

## Script index

| Script | Purpose |
| --- | --- |
| [alignment/starsolo.sh](alignment/starsolo.sh) | STARsolo alignment and gene counting |
| [snrna/01_import_qc_sct.R](snrna/01_import_qc_sct.R) | Count import, sample metadata, QC and SCTransform |
| [snrna/02_integrate.R](snrna/02_integrate.R) | SCT integration, PCA, UMAP and clustering |
| [snrna/03_cell_cycle.R](snrna/03_cell_cycle.R) | S/G2M scoring and marker audit |
| [snrna/04_classify_ivory.R](snrna/04_classify_ivory.R) | Assay-specific ivory expression grouping |
| [snrna/05_find_markers.R](snrna/05_find_markers.R) | Nucleus-level differential expression with explicit contrast membership |
| [snrna/06_expression_ranking.R](snrna/06_expression_ranking.R) | Aggregate expression and CPM ranking |
| [snrna/07_restore_saved_metadata.R](snrna/07_restore_saved_metadata.R) | Restore metadata from individual saved subsets |
| [plots/lineage_panels.R](plots/lineage_panels.R) | UMAP, expression and cell-proportion panels |
| [plots/volcano.R](plots/volcano.R) | Differential-expression volcano plots |
| [plots/dotplot.R](plots/dotplot.R) | Gene-expression dot plots |
| [plots/go_enrichment.R](plots/go_enrichment.R) | Visualization of exported enrichment results |
| [locus/extract_spliced_reads.sh](locus/extract_spliced_reads.sh) | Strand-specific spliced-read extraction |
| [locus/plot_barcode_reads.R](locus/plot_barcode_reads.R) | Barcode-level locus read and ribbon plots |
| [targets/rnahybrid.sh](targets/rnahybrid.sh) | RNAhybrid prediction launcher |
| [targets/summarize_rnahybrid.py](targets/summarize_rnahybrid.py) | Prediction parsing and gene-level summaries |
| [targets/annotate_targets.py](targets/annotate_targets.py) | Comparative DE/target annotation and conserved candidates |

The [analysis index](../README.md#analysis-workflows) describes inputs, settings and paper connections for each workflow. Original source files are retained in `archive/` for provenance.

## Inputs and outputs

Use full Seurat cell names, including sample prefixes, consistently across metadata and membership files. Sample TSVs require `sample_id`, `counts_dir`, `time`, `race` and `genotype`; a header template is in [samples.template.tsv](../config/samples.template.tsv). Contrast TSVs have `barcode` and `group` columns. Gene-label TSVs have `gene_id` and `label` columns. Barcode selection files contain one full cell name per line without a header.

RDS and `.h5seurat` object inputs are supported. R scripts export settings and `sessionInfo()` alongside results. Integration requires an explicit seed. Required numerical settings should come from the matching analysis branch, not from sample filenames. See [the saved-workspace audit](../docs/SAVED_WORKSPACE_AUDIT.md) for subset-specific memberships and assay records.

The differential-expression script implements nucleus-level `FindMarkers`; aggregate expression ranking is a separate analysis. GO plotting consumes exported enrichment statistics and does not calculate a new enrichment test. Multi-layer Seurat objects may require a version-compatible preparation step; layers are not silently joined.

## miR-193 target analysis

The archived predictions can be summarized without rerunning RNAhybrid:

```bash
python3 scripts/targets/summarize_rnahybrid.py \
  --input archive/received/pegasus_miR/Hmel_miR193_RNAhybrid_longest3utr.out \
  --output-dir results/miR193_summary

python3 scripts/targets/annotate_targets.py \
  --annotations data/paper/Table_S8.tsv \
  --summary results/miR193_summary/gene_hit_summary.tsv \
  --output-dir results/miR193_comparative
```

See [the target workflow](../analyses/10_mir193_targets.md) for the prediction command, identifier rules, energy-cutoff distinction and validation against Table S8.

## Software

- Bash, awk, STAR and samtools for alignment and BAM processing.
- R with Seurat, sctransform and SeuratDisk; edgeR for CPM; ggplot2 and ggrepel for expression plots.
- Rsamtools, GenomicRanges, IRanges and scico for the locus plotter.
- Python 3.9 or later for target summarization and annotation.
- RNAhybrid in a micromamba environment for new target predictions.

Exact historical software versions require the matching environment records. Checks against saved tables validate those data comparisons; they do not establish that every scientific workflow has been rerun. See [VALIDATION.md](../VALIDATION.md).
