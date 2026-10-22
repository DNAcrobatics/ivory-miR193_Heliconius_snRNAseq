#!/usr/bin/env Rscript
# Export existing annotations; does not recluster or change the Seurat object.
x <- commandArgs(trailingOnly=TRUE)
if (length(x) != 6L) stop(paste(
  'Usage: Rscript export_seurat_metadata.R object.rds cells.tsv',
  'sample_column celltype_column barcode_column_or_rownames exact_or_last-underscore'))
if (file.exists(x[2])) stop('Output exists')
suppressPackageStartupMessages(library(SeuratObject))
obj <- readRDS(x[1])
meta <- obj[[]]
if (!all(x[3:4] %in% names(meta))) stop('Sample/cell-type columns missing')
bc <- if (x[5] == 'rownames') rownames(meta) else {
  if (!x[5] %in% names(meta)) stop('Barcode column missing')
  as.character(meta[[x[5]]])
}
if (!x[6] %in% c('exact','last-underscore')) stop('Unknown barcode mode')
if (x[6] == 'last-underscore') bc <- sub('.*_', '', bc)
out <- data.frame(sample=as.character(meta[[x[3]]]), barcode=bc,
                  celltype=as.character(meta[[x[4]]]))
if (anyNA(out) || any(out == '') || anyDuplicated(out[c('sample','barcode')]))
  stop('Missing values or duplicate sample/barcode pairs')
dir.create(dirname(x[2]), recursive=TRUE, showWarnings=FALSE)
write.table(out, x[2], sep='\t', quote=FALSE, row.names=FALSE)
cat(nrow(out), 'cells exported to', x[2], '\n')
