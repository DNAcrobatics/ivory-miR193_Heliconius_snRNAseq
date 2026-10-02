#!/usr/bin/env Rscript
# Clean the complete race-specific RNA classification branch.
# Requires explicit thresholds; incomplete RNA_by_Ivory experiment is not guessed.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages(library(Seurat))
a <- parse_args(c("input", "cells", "assay", "iv-threshold", "other-threshold", "output"))
prepare_output(a$output); obj <- read_object(a$input)
# These are exact full Seurat cell names, not raw barcodes without sample prefixes.
cells <- readLines(a$cells); cells <- cells[nzchar(cells)]
if (!length(cells) || anyDuplicated(cells) || any(!cells %in% colnames(obj))) stop("Invalid contrast barcode selection")
obj <- subset(obj, cells = cells)
if (a$assay != "RNA") stop("This recovered branch uses RNA; other assays require separate confirmed workflow")
DefaultAssay(obj) <- "RNA"; check_columns(obj[[]], "race")
if (anyNA(obj$race) || any(!nzchar(obj$race))) stop("Missing race labels")
iv <- FetchData(obj, vars = "Ivory", slot = "data")[colnames(obj), "Ivory"]
if (any(!is.finite(iv))) stop("Missing or non-finite Ivory expression")
threshold <- ifelse(obj$race == "iv", number(a[["iv-threshold"]], "iv-threshold"), number(a[["other-threshold"]], "other-threshold"))
obj$Iv <- ifelse(iv > threshold, "hi", "lo")
obj$ivory_classification_value <- iv; obj$ivory_classification_threshold <- threshold
write_tsv(cbind(barcode = colnames(obj), obj[[]]), paste0(a$output, ".metadata.tsv"))
save_run(obj, a$output, a)
