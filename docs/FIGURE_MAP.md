# Figure to source map after intake

Scope follows `FinalEdits_Ivory_d78k_scRNAseq_.docx` and `Supplementary_dta.docx` supplied 2026-10-01. The larger pasted R version is retained as `archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r`. It contains candidate branches for the following figures; final output matching is pending.

| Figure/panels | Recovered candidate source | Still needed |
|---|---|---|
| 1–4, S1–S5 | Manuscript/supplement legends; image transformation descriptions | Manual assembly/provenance; image code/macros only if actually used |
| 5A–D | SOP subsetting, integration candidates, UMAPs, ivory features and RaceClusters | Combined-object builder, final barcode/annotation table and confirmed plotting branch |
| 5E–F | Early/late subset, high/low classification candidates, PhaseClustersEarly/Late | Final assay/threshold, cell metadata and exported panels |
| 6A | Partial comments/commands in `ivory_snuc/samtools.sh` | Barcode-aware strand/depth pipeline, locus models, coverage tables and ribbon plotter |
| 6B–C | Manuscript transcript/schematic descriptions | Coordinates and manual assembly provenance |
| 6D–E | Expanded paste: `HiVsLo.markers`, early volcano/dotplot code | Saved full DE result, final group memberships and inputs |
| 6F–G | Expanded paste: `HiVsLoLate.markers`, late volcano/dotplot code | Saved full DE result, final group memberships and inputs |
| 7A–B | Expanded paste: `WtVMutInHi.Markers` / `Markers2`, early plot and miRNA highlight code | Resolve variable version, target-hit derivation and full DE table |
| 7C | Expanded paste: imported Up/Down PANGEA DirectGO CSVs, `GoTermsOut.pdf` | Actual CSVs, generating job/settings and correction choice |
| 7D–E | Expanded paste: `WtVMutInHiLate.Markers`, late plot/miRNA highlight code | Full DE table, hit lists and final memberships |
| 7F | Expanded paste: imported late PANGEA DirectGO CSVs, `GoTermsOutLate.pdf` | Actual CSVs, generating job/settings and correction choice |
| S6 | Nanopore Methods only | Dorado/mapping/depth/plot code and inputs |
| S7 | Expanded paste: whole-object marker/dimension/feature plots | Final full-object builder and annotations |
| S8 | Expanded paste: AggregateExpression, edgeR CPM, ranking and histogram | Final source object/cell set, output and binning confirmation |
| Tables S4–S8 and target comparison | Manuscript describes DE, overlap and RNAhybrid results | Machine-readable tables, RNAhybrid/overlap code and external inputs |

Find literal input/export filename references in `metadata/literal_file_references.tsv`. `metadata/figure_manifest.tsv` records candidates, not executable reproduction commands. RNAhybrid/cross-species target analysis is now a required branch because it appears explicitly in the supplied Results/Methods and Figure 7 legend.
