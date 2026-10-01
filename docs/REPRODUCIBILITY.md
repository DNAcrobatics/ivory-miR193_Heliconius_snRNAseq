# Assembly and verification plan

## Preserve the scientific record

Collect original source first. Record original path, content checksum and how each file maps to the final repository. Keep uncertain/obsolete versions in the intake bundle while selecting the version that produced each accepted result. Do not treat a recent filename as proof that it is final. Reorganization may change paths and add configuration, but must not silently change filtering, model definitions, normalization, tests or plotting inputs.

## Run order to document

1. Resolve references, annotation edits, identifiers and orthology.
2. Process snRNA-seq reads with the actual STARsolo configuration and import counts.
3. Apply QC, integration/normalization and annotation in the order established by final code; select SOP-lineage cells and compute cell-cycle summaries.
4. Define early/late membership and high/low ivory groups; execute each exact DE contrast.
5. Run the retained enrichment workflow on its explicit gene universe.
6. In a separate dependent branch, select relevant cell barcodes from alignments and produce strand-specific locus coverage.
7. Run Nanopore genotyping as its own branch.
8. Run the target prediction/cross-species comparisons described in the supplied manuscript.
9. Generate figure panels and tables, then record manual assembly steps.

This is a dependency plan, not an executable workflow. Exact commands and entry points are pending source inspection.

## Two reproducibility levels

- **Plot reproduction:** frozen per-panel tables and metadata, a pinned plotting environment, explicit commands and expected panel outputs.
- **Analysis reproduction:** raw data accessions, reference checksums, exact upstream commands and environments, sample metadata, scripts and validation summaries.

Large raw inputs and computed indexes stay outside Git. The repository records where to obtain them. Large intermediate objects are only required when a documented reproduction entry point depends on them.

## Release criteria

- Every main/supplementary analysis panel and derived table has a manifest entry; manual image panels are explicitly marked.
- Every retained script has an owner/source and input/output description, and resolves project paths via documented configuration or arguments.
- DE replication unit, test, contrasts and multiple-testing procedures are established from code, not inferred from old prose.
- Enrichment implementation, universe, mappings and parameters match the accepted results.
- Seeds, package/tool versions and reference/annotation versions are recorded.
- All branches receive relevant syntax/input checks; at least plot reproduction is run on the supplied inputs and compared to the final outputs. Full reruns are reported separately with their actual validation status.
- Citation, license choice, code release version and archival DOI are completed before the paper's code-availability statement is finalized.

## Current intake status

Original snapshots are in `archive/`. The main DE candidates operate on nucleus subsets; expression aggregation blocks belong to S8. Do not invent a pseudobulk DE branch to reconcile prose. Resolve the final branch, statistical description and group membership using saved results and per-cell metadata. PANGEA DirectGO CSV plotting has been recovered, but enrichment generation/settings are missing. See `INTAKE_REVIEW.md` for the specific blockers.
