#!/usr/bin/env Rscript
# Clean cell-cycle scoring from ivory_reanalysis.expanded.pasted.r.
# The source table is unchanged. Effective marker lists and exclusions are exported.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages(library(Seurat))
a <- parse_args(c("input", "genes", "assay", "output")); prepare_output(a$output)
obj <- read_object(a$input)
if (!a$assay %in% Assays(obj)) stop("Requested assay missing")
DefaultAssay(obj) <- a$assay
mapping <- read_tsv(a$genes); check_columns(mapping, c("stage", "HmID"))
mapping$matched <- !is.na(mapping$HmID) & mapping$HmID %in% rownames(obj)
mapping$repeated_within_phase <- duplicated(mapping[c("stage", "HmID")])
write_tsv(mapping, paste0(a$output, ".marker_audit.tsv"))
# Make absent features and duplicate treatment explicit instead of relying on warnings.
features <- function(stage) unique(mapping$HmID[mapping$stage == stage & mapping$matched])
s <- features("S"); g <- features("G2/M")
if (!length(s) || !length(g)) stop("No usable markers for one or both phases")
obj <- CellCycleScoring(obj, s.features = s, g2m.features = g, set.ident = TRUE)
write_tsv(cbind(barcode = colnames(obj), obj[[]]), paste0(a$output, ".metadata.tsv"))
save_run(obj, a$output, a)
