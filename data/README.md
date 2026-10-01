# Data layout

`processed/` may contain small frozen plotting inputs; `external/` is for locally downloaded external inputs and is ignored by Git. Populate `metadata/data_manifest.tsv` with actual public accessions/URLs, versions, hashes and intended local paths. Input filenames listed in manifests do not imply that the files are included.

Raw reads, BAM/CRAM files, Seurat objects and indexes should be obtained from documented archives or supplemental storage. The first intake pass needs only their manifest and the immediate small tables used by the figures.
