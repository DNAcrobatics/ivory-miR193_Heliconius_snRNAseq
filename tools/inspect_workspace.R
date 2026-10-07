#!/usr/bin/env Rscript
# Inspect saved R objects without sourcing historical analysis or startup code.
# Run with Rscript --vanilla. Large Seurat workspaces may require substantial RAM.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript --vanilla tools/inspect_workspace.R input.RData output_directory")
if (!file.exists(args[1])) stop("Input does not exist")
if (file.exists(args[2])) stop("Output directory already exists")
dir.create(args[2], recursive = TRUE)
# Loading class definitions makes saved Seurat objects inspectable when installed.
if (requireNamespace("Seurat", quietly = TRUE)) loadNamespace("Seurat")
env <- new.env(parent = emptyenv())
if (grepl("\\.rds$", args[1], ignore.case = TRUE)) {
  assign("saved_object", readRDS(args[1]), envir = env)
} else {
  load(args[1], envir = env)
}
object_names <- ls(env, all.names = TRUE)
inventory <- lapply(seq_along(object_names), function(i) {
  n <- object_names[i]; x <- get(n, envir = env, inherits = FALSE)
  dims <- tryCatch(dim(x), error = function(e) NULL)
  # Export ordinary tables (including FindMarkers results), never evaluate
  # stored functions or recursively traverse an arbitrary environment.
  if (is.data.frame(x)) write.table(x, file.path(args[2], paste0("object_", i, ".tsv")), sep = "\t", quote = TRUE, row.names = TRUE)
  if (isS4(x) && "meta.data" %in% slotNames(x)) {
    md <- slot(x, "meta.data")
    if (is.data.frame(md)) write.table(md, file.path(args[2], paste0("object_", i, "_metadata.tsv")), sep = "\t", quote = TRUE, row.names = TRUE)
  }
  data.frame(index = i, object = n, class = paste(class(x), collapse = ";"),
             dimensions = paste(dims, collapse = "x"), bytes = as.numeric(object.size(x)))
})
if (length(inventory)) write.table(do.call(rbind, inventory), file.path(args[2], "objects.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[2], "inspection_sessionInfo.txt"))
writeLines(c("Input:", normalizePath(args[1]), "MD5:", unname(tools::md5sum(args[1])),
             "This session describes the inspection environment, not the original analysis environment."), file.path(args[2], "provenance.txt"))
