# Shared command-line and file utilities. Uses base R; packages load in each script.
parse_args <- function(required) {
  x <- commandArgs(trailingOnly = TRUE)
  if (!length(x) || length(x) %% 2L != 0L) stop("Supply --name value pairs; see scripts/README.md")
  keys <- x[seq(1L, length(x), 2L)]
  if (any(!startsWith(keys, "--")) || anyDuplicated(keys)) stop("Invalid or duplicate option")
  out <- as.list(x[seq(2L, length(x), 2L)]); names(out) <- substring(keys, 3L)
  missing <- setdiff(required, names(out))
  if (length(missing)) stop("Missing options: ", paste(missing, collapse = ", "))
  out
}
number <- function(x, name, minimum = 0) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || !is.finite(value) || value < minimum) stop("Invalid ", name)
  value
}
integer_option <- function(x, name, minimum = 1) {
  v <- number(x, name, minimum)
  if (v != floor(v)) stop(name, " must be an integer")
  as.integer(v)
}
check_columns <- function(x, cols) {
  missing <- setdiff(cols, colnames(x))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}
read_tsv <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
prepare_output <- function(path) {
  if (file.exists(path)) stop("Output already exists: ", path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
read_object <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  if (!grepl("\\.h5seurat$", path, ignore.case = TRUE)) stop("Input must be RDS or h5seurat")
  SeuratDisk::LoadH5Seurat(path)
}
save_run <- function(object, path, options) {
  saveRDS(object, path)
  writeLines(capture.output(dput(options)), paste0(path, ".options.R"))
  writeLines(capture.output(sessionInfo()), paste0(path, ".sessionInfo.txt"))
}
