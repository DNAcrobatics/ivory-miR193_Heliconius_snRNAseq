# ivory:mir-193 and butterfly scale differentiation

Analysis code accompanying *The gene ivory:mir-193 controls scale type differentiation in Heliconius butterflies* (Nature Communications).

This repository contains annotated workflows for single-nucleus RNA sequencing, scale-cell classification, differential expression, gene ontology visualization, ivory locus read visualization and miR-193 target prediction. The scripts use explicit input paths and analysis parameters, and export machine-readable results alongside figures.

## Analysis workflows

| Workflow | Description |
| --- | --- |
| [Reference preparation](analyses/01_reference_preparation.md) | Genome and gene annotation resources |
| [Alignment and counting](analyses/02_starsolo_counts.md) | STARsolo single-nucleus gene counts |
| [Quality control and integration](analyses/03_snrna_qc_integration.md) | Import, filtering, SCTransform and integration |
| [Scale-cell lineage](analyses/04_sop_lineage_annotation.md) | Lineage selection and cell-type annotation |
| [Cell cycle](analyses/05_cell_cycle.md) | S and G2/M scoring and phase visualization |
| [Ivory expression groups](analyses/06_ivory_high_low.md) | Definition of ivory-high and ivory-low nuclei |
| [Differential expression](analyses/07_differential_expression.md) | Early and late expression comparisons |
| [Gene ontology](analyses/08_go_enrichment.md) | Gene mapping and enrichment visualization |
| [Ivory locus](analyses/09_locus_readthrough.md) | Strand-specific spliced reads and barcode-level plots |
| [miR-193 targets](analyses/10_mir193_targets.md) | RNAhybrid predictions and comparative target annotations |
| [Expression ranking](analyses/12_expression_ranking.md) | Aggregate expression and ivory ranking |

Nanopore genotyping is described in [the separate workflow overview](analyses/11_nanopore_crispr_genotyping.md). Its processing was performed by a collaborator; this collection does not include an executable reconstruction of that analysis.

## Getting started

1. Read the workflow page for the analysis of interest.
2. Install the required software described in [the script guide](scripts/README.md) and [environment documentation](envs/README.md).
3. Obtain the workflow inputs and select the corresponding sample metadata, assay and cell membership.
4. Run the commands in the script guide from the repository root, writing each analysis branch to a separate output directory.

The single-nucleus workflows use R with Seurat, sctransform and SeuratDisk. Alignment and locus processing use STAR and samtools. The miR-193 summary and annotation scripts use Python 3.9 or later and its standard library; prediction requires RNAhybrid.

## Repository contents

| Directory | Contents |
| --- | --- |
| [`scripts/`](scripts/README.md) | Annotated analysis scripts and execution examples |
| [`analyses/`](analyses/) | Workflow descriptions and paper connections |
| [`data/`](data/) | Small reference inputs and supplementary tables |
| [`metadata/`](metadata/) | Sample information, file manifests and comparison records |
| [`envs/`](envs/README.md) | Software requirements |
| [`archive/`](archive/README.md) | Original source files retained for provenance |
| [`tools/`](tools/) | Utilities for exporting saved analysis objects |

Large sequencing files, BAMs and complete Seurat workspaces are not stored in this repository. Sequencing resources are associated with BioProjects [PRJNA1404329](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1404329) and [PRJNA1118807](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1118807). Input locations and checksums are recorded in [the data manifest](metadata/data_manifest.tsv).

## Reproducibility

Original source files and distinct analysis branches are preserved. The scripts expose parameter choices rather than inferring them from filenames or figure captions. Saved differential-expression results have been compared with supplementary tables, and the RNAhybrid parser reproduces the supplied gene summaries. These data comparisons do not constitute a complete rerun of the R analysis or sequencing pipeline. See [validation](VALIDATION.md) and the individual workflow pages for the scope of each check.

The miR-193 workflow retains the supplied RNAhybrid energy cutoff of −25 kcal/mol. The manuscript describes −20 kcal/mol; these settings must be reconciled before interpreting the archived predictions as a rerun of that Methods description.

## Citation and reuse

Please cite the accompanying paper and identify the repository revision used. Machine-readable citation metadata is available in [CITATION.cff](CITATION.cff). A software license and archived release DOI have not yet been assigned.
