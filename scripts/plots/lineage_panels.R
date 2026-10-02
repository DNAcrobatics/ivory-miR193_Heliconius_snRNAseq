#!/usr/bin/env Rscript
# Clean UMAP, Ivory feature, and cell-type/phase proportion panels.
# Inputs must already carry confirmed celltypes, Phase, race/time/genotype labels.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
a <- parse_args(c("input", "assay", "group-by", "category", "feature", "output-dir", "seed"))
if (dir.exists(a[["output-dir"]])) stop("Output directory already exists")
obj <- read_object(a$input); check_columns(obj[[]], c(a[["group-by"]], a$category))
if (!is.null(a$cells)) {
  cells <- readLines(a$cells); cells <- cells[nzchar(cells)]
  if (!length(cells) || anyDuplicated(cells) || any(!cells %in% colnames(obj))) stop("Invalid barcode selection")
  obj <- subset(obj, cells = cells)
}
if (!a$assay %in% Assays(obj)) stop("Assay missing")
DefaultAssay(obj) <- a$assay
if (!"umap" %in% Reductions(obj) || !a$feature %in% rownames(obj)) stop("UMAP or feature missing")
meta <- obj[[]]; if (anyNA(meta[c(a[["group-by"]], a$category)])) stop("Missing plotted metadata")
set.seed(integer_option(a$seed, "seed", 0))
dir.create(a[["output-dir"]], recursive = TRUE)
output <- function(name) file.path(a[["output-dir"]], name)
p <- DimPlot(obj, group.by = a[["group-by"]], shuffle = TRUE) + theme_void()
ggsave(output("umap.pdf"), p, device = "pdf", width = 5, height = 4)
f <- FeaturePlot(obj, features = a$feature, order = TRUE) + theme_void()
ggsave(output("feature.pdf"), f, device = "pdf", width = 5, height = 4)
# Denominator is ALL selected nuclei in each group, not only one cell type.
tab <- as.data.frame(table(group = meta[[a[["group-by"]]]], category = meta[[a$category]]))
tab$fraction <- tab$Freq / ave(tab$Freq, tab$group, FUN = sum)
p <- ggplot(tab, aes(x = 1, y = fraction, fill = category)) +
  geom_col() + coord_polar("y") + facet_wrap(~group) + theme_void()
ggsave(output("proportions.pdf"), p, device = "pdf", width = 7, height = 3)
write_tsv(tab, output("proportions.tsv"))
write_tsv(cbind(barcode = colnames(obj), meta), output("cell_metadata.tsv"))
write_tsv(cbind(barcode = rownames(Embeddings(obj, "umap")), Embeddings(obj, "umap")), output("umap_coordinates.tsv"))
writeLines(capture.output(dput(a)), output("options.R"))
writeLines(capture.output(sessionInfo()), output("sessionInfo.txt"))
