DATA AND INPUT GUIDE — ivory:mir-193 paper

This repository contains code, small supplied text inputs and provenance records.
It does not contain raw reads, BAM files, full Seurat objects or private manuscripts.

Public sequencing BioProjects cited during intake:
  PRJNA1404329 — https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1404329
  PRJNA1118807 — https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1118807
Exact run-to-sample associations, download paths and checksums remain to recover.

Small supplied input:
  archive/original/ivory_proj_transfer_aug2025/CellCycleGenes.txt
  Columns: stage, fbgn, code, CeuID, HmID.
  Native identifiers, duplicate mappings and NO_HIT values are preserved.

Input inventory and metadata:
  metadata/data_manifest.tsv — expected inputs and acquisition status.
  metadata/literal_file_references.tsv — filenames recovered from source.
  metadata/sample_metadata_legacy_May2024.tsv — historical values only.
  metadata/sample_manifest.tsv — final metadata template, pending confirmation.

Required external inputs include the edited reference/GTF, filtered count matrices,
final integrated/lineage objects, barcode-tagged alignment files, enrichment exports,
miR-193/UTR sequences and predictions, and Nanopore genotyping inputs/results.
See each numbered workflow for its prerequisites and docs/FILE_REQUESTS.md for the
complete request list. Large files should be deposited externally with stable
locations and checksums; their absence is not represented as a completed archive.

Sample naming and developmental timing:
Do not infer final genotype, race or stage from filenames. Historical time values
differ from revised paper values. Confirm a final sample-level and per-cell mapping.

Coordinate and annotation conventions:
Record genome assembly, scaffold, strand, interval convention and annotation version
for every locus plot and extracted sequence. BED and samtools intervals have different
coordinate conventions; the final conversion rules must be documented with the code.
