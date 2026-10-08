#!/usr/bin/env Rscript
# Restore a frozen saved cohort and its labels without recomputing thresholds.
# This is an audit/recovery step, not proof that a saved DE table used this state.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages(library(Seurat))
a <- parse_args(c("input", "metadata", "assay", "identity", "output"))
prepare_output(a$output)
obj <- read_object(a$input)
meta <- read_tsv(a$metadata)
check_columns(meta, c(".saved_rowname", a$identity))
cells <- meta$.saved_rowname
if (!length(cells) || anyNA(cells) || anyDuplicated(cells) || any(!cells %in% colnames(obj))) {
  stop("Saved cells must be unique and present in the input object")
}
if (!a$assay %in% Assays(obj)) stop("Requested assay is absent; metadata cannot reconstruct expression matrices")
if (anyNA(meta[[a$identity]])) stop("Saved identity contains missing values")
obj <- subset(obj, cells = cells)
rownames(meta) <- cells
meta$.saved_rowname <- NULL
# Retain QC/count fields computed from the input; restore biological assignments.
# Fail if the caller supplies a different library identity for the same cell.
if ("orig.ident" %in% names(meta) && any(as.character(obj$orig.ident) != as.character(meta[colnames(obj), "orig.ident"]))) {
  stop("Input library identities disagree with the saved metadata")
}
fields <- intersect(c("time", "race", "sop", "celltypes", "Iv", "mut", "Phase",
                      "S.Score", "G2M.Score", "old.ident", "norm_factor", ".saved_active_ident"), names(meta))
for (field in fields) obj[[field]] <- meta[colnames(obj), field, drop = FALSE]
DefaultAssay(obj) <- a$assay
Idents(obj) <- setNames(as.character(meta[[a$identity]]), cells)[colnames(obj)]
write_tsv(cbind(barcode = colnames(obj), obj[[]]), paste0(a$output, ".metadata.tsv"))
save_run(obj, a$output, a)
