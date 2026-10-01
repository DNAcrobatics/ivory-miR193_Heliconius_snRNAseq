# ivory:mir-193 — paper analysis code

Code collection for **The gene ivory:mir-193 controls scale type differentiation in Heliconius butterflies** (Hanly and colleagues; Nature Communications; final citation pending).

**Status: original source collected in part; publication release not yet runnable.** Five source snapshots and the cell-cycle input table are preserved in `archive/`. Source provenance, unresolved dependencies and manuscript/code differences are documented. No biological result has been reproduced yet. See [the intake review](docs/INTAKE_REVIEW.md).

## Organization

| Location | Purpose |
|---|---|
| `analysis/01_references/` | Reference preparation, annotation edits, orthology and gene identifiers |
| `analysis/02_snrna/` | STARsolo, QC, integration, clustering, annotation and cell-cycle scoring |
| `analysis/03_differential_expression/` | Exact early/late and high/low ivory comparisons and tests |
| `analysis/04_enrichment/` | Final GO/phenotype enrichment workflow and plotting |
| `analysis/05_locus_readthrough/` | Barcode selection, strand-specific coverage and locus tracks |
| `analysis/06_candidate_targets/` | miRNA target prediction and cross-species overlaps, if used in final paper |
| `analysis/07_crispr_genotyping/` | Nanopore basecalling, mapping, coverage and genotype summaries |
| `analysis/08_figures/` | Figure exporters and any programmatic final assembly |
| `metadata/` | Sample, contrast, figure and provenance manifests |
| `config/` | Portable path configuration templates |
| `envs/` | Actual environment exports and software version records, once supplied |
| `data/` | Small reproducibility inputs, with accession manifests for external data |
| `docs/` | Collection checklist, figure map, execution plan and release criteria |
| `tools/` | Utilities for assembling this archive |
| `archive/` | Received originals and exploratory versions, preserved with provenance |

Start with [the file request checklist](docs/FILE_REQUESTS.md). [The figure map](docs/FIGURE_MAP.md) links the earlier manuscript panels to the expected analysis components. It now follows the supplied `FinalEdits_Ivory_d78k_scRNAseq_.docx` and supplementary figures; final panel/output matching remains pending.

## Reproduction plan

We will support two entry points: (1) regenerate paper plots and tables from frozen processed inputs, and (2) rerun upstream processing from public raw reads and pinned references. Each component will declare inputs, outputs, software versions, parameters and run order. Cluster submission scripts will accompany the underlying commands; readers should not need your HPC account to understand the computation.

The intended scientific order is documented in [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md). Runnable commands will be added only after the original source and its inputs have been inspected. Existing numerical script names and originals will retain a provenance mapping when reorganized.

## Data

The supplied manuscript and reporting summaries list snRNA-seq BioProjects **PRJNA1404329** and **PRJNA1118807**. These are manuscript-reported identifiers, not yet verified download manifests. Individual run accessions and their sample mappings are still needed. Large reads, BAMs, reference indexes and analysis objects will be stored externally with documented accessions or download locations and checksums. Small required lookup tables and plotting inputs can live here.

## Citation and license

`CITATION.cff` contains a provisional paper title and author entry. Replace with the final citation and code-release identifier before publication. A code license has not been selected; see `docs/LICENSE_DECISION.md`. Publication of this scaffold does not itself grant an open-source license.

## Repository

https://github.com/Hanliconius/ivory-mir193_paper

## Latest intake

The next requests are the producer of `Hm_combined.sct.h5seurat`, the exact small lookup/GO files named in the recovered code, final sample/per-cell metadata and saved DE/source-data tables. Read [FILE_REQUESTS.md](docs/FILE_REQUESTS.md) for the precise filenames and a batch collection command. The original R snapshots contain incomplete/session-dependent code and are not offered as release entry points. The source manuscripts were inspected but are not included in this code archive.
