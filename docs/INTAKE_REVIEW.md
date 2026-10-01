# Intake review on 2026-10-01

The repository now contains five source snapshots and one cell-cycle input table. The two attached text originals are retained unchanged. Manuscript/editorial documents were inspected as scope evidence; they are not copied into this code archive.

## What the source establishes

The longer paste has candidate code for Figure 5, Figure 6D–G, Figure 7, S7 and S8. Main early scales use `iv_20`, `iv_61`, `Hmm_25`, `Hmr_25`, `Hmr_30`; the later main subset uses `iv_55`, `Hmm_40`, `Hmr_40`. Candidate DE calls use `FindMarkers` on cell subsets with `test.use` omitted, and plot `avg_log2FC` and `p_val_adj`. The supplied legends specify Wilcoxon rank-sum/Bonferroni. Exploratory `test.use='roc'` calls elsewhere are separately retained and are not assumed to generate the published p values.

The supplied `AggregateExpression`/edgeR CPM blocks rank expression for the supplementary histogram. They do not establish a pseudobulk DE analysis. Complete output matching is needed before selecting a release branch.

PANGEA-named CSV files are imported to draw DirectGO-category plots. No actual enrichment-generating job or algorithm is supplied. The permutation claim comes from the manuscript legend, not from the plot code. BY and BH column names appear in differing filter/color choices. Request the source exports/settings rather than inventing a custom enrichment implementation.

RNAhybrid and cross-species target comparisons are now confirmed parts of the manuscript. The source supplies a 15-gene highlight list but not its derivation or the `with_miRhit` definition. The Methods give RNAhybrid 2.1.2, `-s 3utr_fly` and `-e -20`; actual command/input provenance remains missing.

## Main execution gaps

- `Hm_combined.sct.h5seurat` is loaded at the start; its producing script is missing. The May 2024 source saves a differently named `Iv_combined.full.sct` and does not import all time-series samples used later.
- Both R snapshots contain literal prose and incomplete code. The longer one includes an unfinished `FindNeighbors` call, shell/awk commands and an incomplete `Ivory > 0.5 &` expression. It also references session variables before their definitions.
- High/low definitions differ across branches: raw/log-normalized RNA thresholds vs an ivory-normalized assay. The assay-building block appears after several analyses and is followed by incomplete assignment code; final branch order cannot be inferred safely.
- `WtVMutInHi.Markers2` is computed while neighboring plotting code refers to `WtVMutInHi.Markers`. The saved table is needed to resolve this session dependency.
- `starsolo.sh` includes GeneFull, EM multimapper handling, EmptyDrops_CR, R2/R1 input order and 16-base CB/12-base UMI. Its reference/index and `jobs` are missing; BAM/CB output options are only comments.
- `samtools.sh` passes shell parsing but its `list` assignment is commented out. Active commands extract spliced alignments from already prepared strand BAMs; they do not generate the complete Figure 6A ribbon tracks.
- The May 2024 sample stages differ from the manuscript. Its initial object construction requires 2,000 features before a later 200-feature subset, so the later threshold alone is not the effective initial filter. Integration/PC/resolution settings also differ among versions and prose.
- The cell-cycle mapping has 66 rows, 3 Hmel `NO_HIT` rows, and repeated Hmel IDs. Preserve these originals; only the actual final CellCycleScoring call/warnings can establish the effective feature list.

## Manuscript/release alignment

Results/Methods still call DE pseudobulk; current legends state nucleus-level Seurat Wilcoxon tests. Figure 7's introductory Type I label conflicts with Results' presumptive Type II description. These are documentation discrepancies to resolve once the final code/output branch is established; the manuscript has not been edited.

The reporting PDFs state STAR 2.7.11b, samtools 1.2.1 and Seurat 5.4.0. These are manuscript-reported values and require historical logs/session info. The dynamic `nr-reporting-summary.pdf` was inspected through its embedded XFA dataset; normal PDF text extraction shows only a viewer placeholder. The other reporting PDF has a flattened readable summary.

The manuscript's Code Availability points to the Hanliconius GitHub profile. The author checklist explicitly requests a DOI-minting archive for custom code. Final delivery should use a specific repository plus tagged archived release; no DOI or repository URL has been invented.

## Validation status

Original attachment bytes/checksums and extracted terminal body verified; inline table structure checked; both shell scripts checked with `bash -n`; archive integrity and internal Markdown links checked. R is unavailable in this workspace, so no R parse/runtime test or biological analysis rerun has been performed. Static defects are documented rather than hidden by replacing the source.
