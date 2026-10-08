# Paper documents recovered from Drive — 8 October 2026

The request for manuscript/supplement/reporting documents and final figure locations is closed. These are accessible in Drive, not absent inputs:

- [2024_Ivory_d78k_scRNAseq](https://drive.google.com/drive/folders/1e_iWDpfH-nYLZICY4Op0LfzHXr_KP6yK): resubmission manuscript, cover letters, reviewer response, supp_tables and figures.
- [Editoral_publication_finals](https://drive.google.com/drive/folders/1BrKB7u3ffX2EbUwIlI3NQf7-aAw2m5R6), inside iv_reanalysis: FinalEdits manuscript, Supplementary_dta, checklist and both reporting PDFs. These document contents were already inspected during the initial intake.
- [ALL_FIGURES_NEW](https://drive.google.com/drive/folders/1IShW0tS_U0F4ebjaJqTaSEk4JL6jGSEW): assembled Illustrator/PDF-compatible candidates, including Figures 5–7 and Figure7 V2. Locations are recovered; the large assemblies have not all been downloaded or visually compared, and version names alone do not prove final publication selection.
- [supp_tables](https://docs.google.com/spreadsheets/d/1-9f8XzSKwGWFgAGpr_4uk7QYfzfEkOs6s0wKyc_besk/edit): exported read-only to data/paper/supp_tables_20261008.xlsx, with all eight tabs also exported to TSV. The spreadsheet export preserves the snapshot; TSV conversion uses the second row as headers, retains data rows and excludes only columns proven entirely blank. Captions are in metadata/paper_supplement_tables.tsv.

S4/S5/S6/S7 contain 99/175/196/102 selected DE rows, respectively. All gene IDs occur in the corresponding saved workspace tables (0067/0073/0221/0222); five statistic columns agree at displayed precision (relative tolerance 0.006). Two S6 p-values displayed as zero are nonzero subnormal values (~1e-309) in the saved table. This is a rounded publication-table match, not an exact byte/numeric match or a scientific rerun. See metadata/paper_saved_de_comparison.tsv. The S7 caption repeats the early 25–32% stage while mapping to Fig.7D; retained unchanged pending editorial reconciliation.

S8 has 733 rows, 111 miR193_RNAhybrid_hit=True rows and 15 candidate_conserved_miR193_target=True rows. It supplies the downstream target/overlap annotations and candidate set previously requested. It does not recover the RNAhybrid command, source UTRs or prediction-to-gene preparation/parsing code. The TSV preserves original column names and annotations; do not infer the exact historical with_miRhit object construction from it alone.

The user states Nanopore processing was performed by a collaborator. It is outside the user's code-intake request: document the contribution and source when available, without asking the user to reconstruct that pipeline or claiming it is reproduced.

Remaining paper-related work is matching existing documents/assemblies to code branches, resolving inconsistent Methods/legend descriptions, and choosing code license/contributor credits/release archive. No repeat paper-document or supplementary-table upload is needed.
