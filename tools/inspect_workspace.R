#!/usr/bin/env Rscript
# Fresh-process wrapper: loads the full workspace once. For an already loaded
# R session, source export_workspace_audit.R instead to avoid another full load.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript --vanilla tools/inspect_workspace.R input.RData output_directory")
if (!file.exists(args[1])) stop("Input does not exist")
if (file.exists(args[2])) stop("Output directory already exists")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "export_workspace_audit.R"))
if (requireNamespace("Seurat", quietly = TRUE)) loadNamespace("Seurat")
env <- new.env(parent = emptyenv())
if (grepl("\\.rds$", args[1], ignore.case = TRUE)) {
  assign("saved_object", readRDS(args[1]), envir = env)
} else {
  load(args[1], envir = env)
}
export_workspace_audit(outdir = args[2], envir = env)
writeLines(c("Input:", normalizePath(args[1]), "MD5:", unname(tools::md5sum(args[1])),
             "This session describes inspection, not necessarily the original analysis environment."), file.path(args[2], "input_provenance.txt"))
