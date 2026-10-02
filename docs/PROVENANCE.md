# Repository provenance

Scaffold created 2026-10-01 and updated after source intake the same day for Joseph Hanly's Nature Communications ivory:mir-193 paper.

Scope evidence now comes from `FinalEdits_Ivory_d78k_scRNAseq_.docx`, `Supplementary_dta.docx`, `Author_Checklist_NCOMMS-26-020429A_1789986563_55.docx` and the two supplied reporting PDFs. These documents supersede the earlier manuscript as scope evidence; internal discrepancies are recorded in `INTAKE_REVIEW.md`. They were read without alteration and are not included in the code package.

Two attached text files were copied byte-for-byte, with both versions retained. The terminal transcript body was extracted after removing boundary shell prompts. Three inline source scripts and CellCycleGenes.txt were transcribed; their exact original cluster bytes are not yet available. See `metadata/received_source_manifest.tsv` for hashes and transformations. No scientific logic was changed. Newly authored documentation, metadata manifests and collection utility are distinguished from the original analysis snapshots.

RNAhybrid/cross-species comparisons are confirmed in the supplied manuscript. Reported software versions are labeled as reported, not verified runtime versions. Sample metadata from the historical May 2024 source remains in a separate legacy manifest and has not been promoted to the final sample sheet.

No analysis has been rerun and no public archive release has been created. The collection is being uploaded to the private repository `Hanliconius/ivory-mir193_paper` under the user's instruction. R is unavailable in this execution workspace. Static review and shell/collection/archive checks are the extent of validation so far.

Received `iv_v_wt_timeseries_March2025/rscript.regressingTime.R` as an inline paste and `Iv_v_Wt_comp/wt_vs_ivory_scales.r` as an attached terminal transcript. The latter attachment is retained byte-for-byte; only its first listing and cat prompt are removed in the extracted .r. Internal console prompts, printed results and exploratory branches remain unchanged.

The March source is transcribed with no logic changes. Both new scientific sources were inspected statically; no analysis was run. The file requests now distinguish a recovered integration/save tail from missing upstream setup.

## Cleaned working scripts — 2026-10-02

The main entry points are now [scripts/README.md](../scripts/README.md). Coherent import/QC/SCT, integration, scoring, classification, DE and plot steps have been extracted and annotated, with explicit arguments and documented reconstruction choices. Original snapshots remain unchanged. Remaining requests concern inputs, final branch selection, missing analysis branches and result matching—not a request to send a pre-cleaned version of every existing script. See [CLEANUP_NOTES.md](CLEANUP_NOTES.md).
