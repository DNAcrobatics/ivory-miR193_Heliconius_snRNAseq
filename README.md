# Code for *The gene ivory:mir-193 controls scale type differentiation in Heliconius butterflies*

Companion code collection for the Nature Communications paper. **Start with the analysis index below.** The repository follows the numbered, annotated workflow format of [Luca’s Anticarsia melanism repository](https://github.com/MilesLuca/Anticarsia_melanism).

Each workflow documents inputs, outputs, software, execution and recovered code. **This is an initial collation, not yet a complete reproduction package.** Supplied scripts and session fragments are preserved; missing sources and unresolved final-analysis choices are explicitly marked.

- **Code and execution guide:** [README_code.md](README_code.md)
- **Data and naming conventions:** [README_DATASET.TXT](README_DATASET.TXT)
- **Outstanding files:** [complete file request checklist](docs/FILE_REQUESTS.md)
- **Figure coverage:** [figure-to-code map](docs/FIGURE_MAP.md)
- **Provenance and validation:** [VALIDATION.md](VALIDATION.md), [source checksums](SOURCE_SHA256SUMS.txt)
- **Citation:** [CITATION.cff](CITATION.cff), pending final publication details and archive DOI
- **License:** awaiting your choice; see [license decision](docs/LICENSE_DECISION.md)

## Analysis index

| Annotated workflow | Paper connection | Current availability |
| --- | --- | --- |
| [01. Reference genome and annotation](analyses/01_reference_preparation.md) | Alignment and all gene-based analyses | Source needed |
| [02. Single-nucleus alignment and counting](analyses/02_starsolo_counts.md) | Upstream of Figures 5–7 and S7–S8 | Historical script |
| [03. Count import, QC and integration](analyses/03_snrna_qc_integration.md) | Full dataset underlying Figures 5–7 and S7 | Historical script; final builder needed |
| [04. SOP lineage integration and annotation](analyses/04_sop_lineage_annotation.md) | Figure 5A–D; Supplementary Figure S7 | Session fragments |
| [05. Cell-cycle scoring and phase plots](analyses/05_cell_cycle.md) | Figure 5E/F | Scoring code and gene mapping supplied |
| [06. Define ivory-high and ivory-low scale nuclei](analyses/06_ivory_high_low.md) | Figures 5E/F and 6–7 | Alternative branches; final choice needed |
| [07. Early and late differential expression](analyses/07_differential_expression.md) | Figure 6D–G and Figure 7A/B/D/E | Session fragments; final tests unresolved |
| [08. Gene mapping, GO enrichment and plotting](analyses/08_go_enrichment.md) | Figure 7C/F | Plotting supplied; enrichment provenance needed |
| [09. Strand-specific ivory locus coverage](analyses/09_locus_readthrough.md) | Figure 6A | Partial historical BAM processing |
| [10. Candidate miR-193 target prediction](analyses/10_mir193_targets.md) | Figure 7 target highlighting | Gene list supplied; prediction code needed |
| [11. Nanopore CRISPR genotyping](analyses/11_nanopore_crispr_genotyping.md) | Supplementary Figure S6 | Source needed |
| [12. Aggregate expression and ivory ranking](analyses/12_expression_ranking.md) | Supplementary Figure S8 | Session fragments |

Figures 1–4 and S1–S5 largely involve experiments, images, schematics and manual assembly. Their source material and any computational steps are tracked in the [figure map](docs/FIGURE_MAP.md), rather than represented by invented analysis code.

## Using the code

1. Choose a workflow from the index and read its status and prerequisites.
2. Obtain its inputs and the exact environment. Large reads, BAMs and Seurat objects belong in external data archives, with locations and checksums recorded in [the data manifest](metadata/data_manifest.tsv).
3. Use the linked original script or recovered code block. Blocks currently retain historical paths, variable names and session dependencies; they are not all standalone scripts. Source-needed pages identify the code still to recover.
4. Confirm the final sample metadata, parameter branch and contrast definitions before replacing historical paths or submitting jobs.
5. Save results with their software logs and compare them with the final figure source tables and vector PDFs.

The main dependencies are reference preparation → alignment/counts → full-object integration → lineage annotation → classification/contrasts → plots. Cell-cycle scoring, coverage, target prediction and genotyping have separate prerequisites documented on their pages.

## Data and naming conventions

Public sequencing resources cited during intake: [PRJNA1404329](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1404329) and [PRJNA1118807](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1118807). Run-to-sample mappings and file checksums still need to be recovered.

Historical sample times and labels in [sample_metadata_legacy_May2024.tsv](metadata/sample_metadata_legacy_May2024.tsv) are preserved as historical values. The [final sample manifest](metadata/sample_manifest.tsv) awaits confirmed metadata; revised paper stages must not be inferred from sample names.

Keep native feature identifiers, including `Ivory`, as used in each assay. Gene-symbol joins, coordinate conventions and custom annotation edits must accompany the final workflows.

## Completing this collection

The [file request checklist](docs/FILE_REQUESTS.md) gives filenames, priorities and the decisions each file resolves. The [collection helper](tools/collect_project_code.py) can collect code and small text inputs from your project directories. Original snapshots remain in [archive/](archive/README.md); manifests and intake notes remain in `metadata/` and `docs/`.

Before release, recover the remaining source, resolve the final analysis choices, validate against paper outputs, choose a license, and archive a versioned release with a DOI. No archive DOI or validated final environment is claimed yet.
