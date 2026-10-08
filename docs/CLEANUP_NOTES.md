# Clean-script preparation

Prepared 2026-10-02 in response to the request for cleaned, annotated scripts. The main code now lives in `scripts/`; original scientific files remain unchanged in `archive/`. The 13 entry scripts and shared base-R helper have explicit inputs and documentation in [scripts/README.md](../scripts/README.md).

## Repairs and adaptations

Removed terminal prompts, printed results, bare notes and abandoned trials from the working code by extracting coherent steps into separate scripts. Original snapshots preserve that material. Replaced implicit session variables and working-directory paths with required arguments. Replaced plotting via an interactive device/dev.copy with explicit ggsave vector PDFs. Added input checks, output checks, recorded settings and sessionInfo; integration seeds are explicit new reproducibility settings rather than claimed original seeds.

The March import/QC/SCT reconstruction comes from commented upstream code and requires sample metadata and numerical settings. The active March integration tail is kept separate so it can also accept the actual prepared Hm_merged object. Metadata labels are supplied; no revised times or final sample membership are guessed. The 2023 glmGamPoi/Phase-regression trials are not mixed into that branch.

Fixed lineage FindNeighbors syntax by using a complete shared integration stage. SOP membership must be a supplied frozen barcode selection; cluster-to-cell-type annotation is not inferred. Cell-cycle scoring writes a marker audit and explicitly excludes unmatched features/deduplicates within phase, a documented adaptation that needs historical comparison. The complete RNA ivory-high/low branch is implemented with explicit thresholds; incomplete RNA_by_Ivory classification is not completed speculatively.

DE requires an explicit two-group membership table, assay and feature-filter values. Wilcoxon is explicit to match the recovered default-style calls/current figure legends; it does not implement or validate the conflicting pseudobulk description. Positive fold changes mean group 1 greater than group 2. Outputs avoid the old Markers/Markers2 collision. Expression aggregation is separate from DE.

GO plotting uses original exported results and explicit filter/color columns, preserving the ability to document differing BY/BH choices without silently selecting a final correction. Volcano display floors zero p-values only for plotting; original values remain in the data. Plots are clean equivalents whose colors, labels, dimensions and selected inputs still need matching to final figures.

STAR counting options and read order are unchanged. STAR output prefixes gain a trailing dot separator. Inactive BAM/CB experiments are not enabled. The locus script implements only the active CIGAR-N step on already strand-separated BAMs. No missing reference-editing, RNAhybrid, Nanopore or per-barcode coverage algorithm has been invented.

## Checks performed

- Bash syntax checks passed for both cleaned shell scripts.
- Mock-tool execution checked STAR sample-array selection, R2/R1 order, retained counting arguments, invalid-index rejection and overwrite rejection.
- Mock-tool execution checked that splice extraction retains headers/CIGAR-N reads, keeps strands separate, creates index outputs and rejects overwrite. These tests verify orchestration, not real BAM processing or STAR alignment.
- Existing source hashes were checked against received_source_manifest.tsv; originals remain unchanged.
- New local Markdown links were checked.

R and the scientific packages are unavailable in this workspace. R scripts have been manually reviewed but have not been parsed by R or executed. Run them in the recovered environment with the actual data, then compare against final DE/source tables and vector PDFs. Missing inputs and competing final parameter branches are tracked in [FILE_REQUESTS.md](FILE_REQUESTS.md).


## Workspace audit additions

The updated inspector delegates to `tools/export_workspace_audit.R`, which can run inside an already loaded session and exports compressed saved tables, Seurat metadata/identities, assay/layer summaries, command records, ggplot inputs and small vectors. Large matrices are excluded. This is recovery tooling, not a replacement analysis. It has been manually reviewed but not executed here; see [WORKSPACE_AUDIT.md](WORKSPACE_AUDIT.md). The refreshed history is unchanged. Saved-object validation remains pending; no history command is treated as proof of successful execution.

## Saved-export recovery — 8 October 2026

The user executed the exporter. Selected outputs pass manifest checksum and DE/plot-input comparisons; see [SAVED_WORKSPACE_AUDIT.md](SAVED_WORKSPACE_AUDIT.md). Added 07_restore_saved_metadata.R to restore each saved subset independently. Parent high/low labels disagree with contrast subsets and must not replace them. Original settings remain unchanged. The helper was reviewed but not run here; data checks are not a scientific rerun.

## Reference/source correction — 8 October 2026

Four reference inputs are now received and structurally checked; the GTF is losslessly compressed with original/stored hashes. The supplied ivory_reanalysis(1).r is recognized as lineage/labeling/DE source. Its early WtVMutInHi.Markers call uses all of test, illustrating why the object name alone cannot prove a high-only contrast. See [REFERENCE_INTAKE.md](REFERENCE_INTAKE.md). Original code and annotations remain unchanged.
