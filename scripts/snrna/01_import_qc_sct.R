#!/usr/bin/env Rscript
# Reconstruct the commented upstream steps in rscript.regressingTime.R.
# Samples, times, genotype and count paths come from a supplied TSV, never names.
# QC order is intentional: min.features at creation precedes the later QC subset.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages(library(Seurat))
a <- parse_args(c("samples", "ribo", "mito", "output", "min-cells", "min-features", "qc-features", "max-mt", "variable-features", "regress", "seed"))
prepare_output(a$output)
samples <- read_tsv(a$samples)
check_columns(samples, c("sample_id", "counts_dir", "time", "race", "genotype"))
if (!nrow(samples) || anyNA(samples) || anyDuplicated(samples$sample_id) || any(!nzchar(samples$sample_id))) stop("Incomplete or duplicate sample metadata")
samples$time <- suppressWarnings(as.numeric(samples$time))
if (any(!is.finite(samples$time))) stop("Times must be numeric")
if (any(!dir.exists(samples$counts_dir))) stop("Count matrix directory missing")
seed <- integer_option(a$seed, "seed", 0); set.seed(seed)
objects <- lapply(seq_len(nrow(samples)), function(i) {
  counts <- Read10X(samples$counts_dir[i])
  if (is.list(counts)) stop("Read10X returned multiple modalities; supply gene-expression directory")
  obj <- CreateSeuratObject(counts, project = samples$sample_id[i],
    min.cells = integer_option(a[["min-cells"]], "min-cells"),
    min.features = integer_option(a[["min-features"]], "min-features"))
  obj$time <- samples$time[i]; obj$race <- samples$race[i]; obj$genotype <- samples$genotype[i]
  obj
})
if (length(objects) < 2L) stop("At least two samples are needed for this merged pipeline")
obj <- merge(objects[[1]], y = objects[-1], add.cell.ids = samples$sample_id, project = "Hm_merged_series")
# Preserve source behavior: build QC gene lists from features shared by all objects.
shared <- Reduce(intersect, lapply(objects, rownames))
genes <- function(path) intersect(shared, read.table(path, header = FALSE, stringsAsFactors = FALSE)[[1]])
ribo <- genes(a$ribo); mito <- genes(a$mito)
if (!length(ribo) || !length(mito)) stop("QC gene lists have no matches; check feature identifiers")
obj <- PercentageFeatureSet(obj, features = ribo, col.name = "percent.Rb")
obj <- PercentageFeatureSet(obj, features = mito, col.name = "percent.MT")
qc_features <- number(a[["qc-features"]], "qc-features")
max_mt <- number(a[["max-mt"]], "max-mt")
obj <- subset(obj, subset = nFeature_RNA > qc_features & percent.MT < max_mt)
if (!ncol(obj)) stop("No nuclei survived QC")
if (any(!samples$sample_id %in% obj$orig.ident)) stop("One or more samples lost all nuclei after QC")
regress <- if (a$regress == "none") NULL else strsplit(a$regress, ",", fixed = TRUE)[[1]]
if (length(regress)) check_columns(obj[[]], regress)
obj <- SCTransform(obj, variable.features.n = integer_option(a[["variable-features"]], "variable-features"),
                   vars.to.regress = regress, seed.use = seed, verbose = TRUE)
write_tsv(cbind(barcode = colnames(obj), obj[[]]), paste0(a$output, ".metadata.tsv"))
save_run(obj, a$output, a)
