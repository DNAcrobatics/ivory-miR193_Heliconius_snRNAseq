# Drive recovery — 7 October 2026

The supplied **Joe** folder was inventoried through its `iv_reanalysis` and `Ivory_Snuc` branches, including analysis, filtered-count, locus and environment subfolders. The initial connector listing exposed many small inputs and PDFs but no saved workspace. **Later audit:** the 15.36 GB `.RData` is now listed; see [WORKSPACE_AUDIT.md](WORKSPACE_AUDIT.md) for verified status, history findings and the compact export needed to inspect its contents. Several nested folders returned no files, including filtered matrix and BAM folders. This records what was accessible; it does not establish that the underlying project lacks those files.

## Recovered material

`metadata/drive_inventory.tsv` records file/folder metadata; `metadata/drive_received_files.tsv` records downloaded files and SHA-256 checksums. Original scripts live in `archive/drive/`; recovered input candidates live in `data/drive/`. These filenames are not declarations that a particular version produced the published panel.

- Two complete `ivory_reanalysis` snapshots, a session history, a larger `wt_vs_ivory_scales.r`, `TS_but_just_SOPs.r`, `millefy_sc_vis.r` and an annotation-extension candidate.
- Four requested early/late GO exports plus older late variants, input gene/protein lists, MODPhea exports, homolog and plot-label lookups, and the cell-cycle list.
- `bcdf.all.txt`: per-barcode cell-type/sample/time mapping useful for the locus plots, not a complete final Seurat metadata export.
- Historical `sample_info.tsv` and sequencing-index sheet. The mutant stages differ from subsequent scripts; retain these as evidence rather than silently replacing final metadata.
- Many requested vector panel PDFs are accessible and catalogued by their Drive links. They are comparison candidates, not yet verified final outputs.

An empty progress-report script was excluded. `Colour_analysis_code.r` identifies itself as code for the unrelated sulphur-butterfly colour/QTL paper and was excluded from this repository.

## Cleaned code recovered from this intake

`scripts/locus/plot_barcode_reads.R` replaces the duplicated sample blocks in `millefy_sc_vis.r` with one annotated command. It extracts alignment starts and CB tags, selects sample-specific metadata, keeps barcodes with more than five alignments, and writes read plots, ribbons, underlying tables and a matching audit. Sample-specific joining avoids collisions when a barcode occurs in multiple libraries. The source loop joined bare barcodes against an all-sample table.

The source ribbon combines 100 successive **occupied start positions** and draws 500 bp rectangles, with an imola log10-count scale limited to 2–3. The cleaned script retains that unusual calculation; it does not call it genomic coverage. Vector PDFs replace interactive `dev.copy` and rasterized marks. Supply the extraction interval used for the BAM: source displays span 1.3–1.7 Mb, while earlier BAM extraction commands span 1.4–1.6 Mb. The annotation track remains separate pending the exact locus GFF. This is a reconstructed candidate, not a validated final Figure 6 reproduction.

`tools/inspect_workspace.R` can inventory a saved workspace and export ordinary result tables and Seurat metadata on a machine with R and the necessary class definitions. It loads objects in an isolated environment and does not source the original analysis history. Its `sessionInfo` describes inspection, not the historical analysis run. The workspace is now listed, but cannot be loaded in this runtime; use the compact export route in [WORKSPACE_AUDIT.md](WORKSPACE_AUDIT.md).

## Still needed, in priority order

1. A compact export from the now-listed `.RData`, using `tools/export_workspace_audit.R` in an already loaded session. Its direct link is recovered; no need to resend it. Export object inventory, saved tables and Seurat metadata as described in [WORKSPACE_AUDIT.md](WORKSPACE_AUDIT.md).
2. Final object-building code, executed parameters and corrected sample metadata. Confirm mutant ages and the integration branch; the recovered historical sample sheet adds another conflicting stage assignment.
3. Final high/low ivory labels, assay and four contrast memberships/full result tables; confirm nucleus-level versus pseudobulk DE. GO CSVs have been recovered, but enrichment settings, tested background and final version selection remain missing.
4. Indexed locus BAMs with CB tags and `Hmel_cortexIvory_region_luca-reannotated.gff`; confirm the final barcode table and extraction interval. The read/ribbon plotting code has now been recovered.
5. RNAhybrid command, sequences/UTR preparation, `with_miRhit` and full targets; Nanopore processing/genotyping code and depth tables; machine-readable supplementary tables.
6. Actual run logs/package versions, final panel-to-code mapping, license and archival release details.

You do not need to resend the recovered GO exports, gene lists, homolog lookup, label tables or historical sample/index sheets. Larger reference annotations and figures can remain externally hosted with checksums and links once their final identity is established.

## Validation

Downloaded file sizes were checked against Drive metadata and each copied file was hashed. CSV headers and row counts were inspected; the early PANGEA CSV contains the columns expected by the existing cleaned GO plotter. R/Seurat execution and scientific result matching remain pending because this environment has no R runtime or count matrices/BAMs, and cannot load the newly listed large workspace. Originals were not run.

