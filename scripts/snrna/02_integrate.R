#!/usr/bin/env Rscript
# Clean March integration/save tail; also supports the supplied SOP integration.
# --cells optionally selects a frozen one-column barcode file for lineage analysis.
# It does not infer SOP clusters or rewrite existing cell-type labels.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages(library(Seurat))
a <- parse_args(c("input", "output", "split-by", "features", "pcs", "dims", "resolution", "seed"))
prepare_output(a$output); obj <- read_object(a$input)
if (!is.null(a$cells)) {
  cells <- readLines(a$cells); cells <- cells[nzchar(cells)]
  if (!length(cells) || anyDuplicated(cells) || any(!cells %in% colnames(obj))) stop("Invalid lineage barcode selection")
  obj <- subset(obj, cells = cells)
}
check_columns(obj[[]], a[["split-by"]])
if (anyNA(obj[[]][[a[["split-by"]]]])) stop("Missing integration group")
if (!"SCT" %in% Assays(obj)) stop("Input lacks SCT assay; run 01_import_qc_sct.R or supply prepared object")
seed <- integer_option(a$seed, "seed", 0); set.seed(seed)
pcs <- integer_option(a$pcs, "pcs"); dims <- integer_option(a$dims, "dims")
if (dims > pcs) stop("dims exceeds computed PCs")
DefaultAssay(obj) <- "SCT"
parts <- SplitObject(obj, split.by = a[["split-by"]])
if (length(parts) < 2L) stop("At least two integration groups required")
features <- SelectIntegrationFeatures(parts, nfeatures = integer_option(a$features, "features"))
parts <- PrepSCTIntegration(parts, anchor.features = features)
anchors <- FindIntegrationAnchors(parts, normalization.method = "SCT", anchor.features = features)
obj <- IntegrateData(anchors, normalization.method = "SCT")
DefaultAssay(obj) <- "integrated"
obj <- RunPCA(obj, npcs = pcs, seed.use = seed)
obj <- RunUMAP(obj, dims = seq_len(dims), seed.use = seed)
obj <- FindNeighbors(obj, dims = seq_len(dims))
obj <- FindClusters(obj, resolution = number(a$resolution, "resolution"), random.seed = seed)
obj <- PrepSCTFindMarkers(obj)
write_tsv(cbind(barcode = colnames(obj), obj[[]]), paste0(a$output, ".metadata.tsv"))
save_run(obj, a$output, a)
# RDS is the portable default. Optional H5Seurat export maintains the historical format.
if (!is.null(a$h5seurat)) {
  prepare_output(a$h5seurat)
  SeuratDisk::SaveH5Seurat(obj, filename = a$h5seurat, overwrite = FALSE)
}
