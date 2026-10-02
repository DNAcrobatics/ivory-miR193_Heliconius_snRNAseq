#!/usr/bin/env Rscript
# Clean DirectGO/PANGEA CSV plotting. Does not generate enrichment statistics.
# Explicit filter/color columns retain each selected historical correction choice.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages(library(ggplot2))
a <- parse_args(c("input", "filter-column", "color-column", "cutoff", "title", "output"))
prepare_output(a$output)
# read.csv deliberately retains the same make.names column conversion as source.
x <- read.csv(a$input, stringsAsFactors = FALSE)
check_columns(x, c("Count.Overlap.Gene", "Gene.Set.Name", "Fold.Enrichment", a[["filter-column"]], a[["color-column"]]))
filter_column <- a[["filter-column"]]; color_column <- a[["color-column"]]
cutoff <- number(a$cutoff, "cutoff")
for (col in c("Count.Overlap.Gene", "Fold.Enrichment", filter_column, color_column)) {
  if (!is.numeric(x[[col]])) stop("Non-numeric result column: ", col)
}
keep <- is.finite(x[[filter_column]]) & x[[filter_column]] < cutoff
x <- x[keep, , drop = FALSE]
if (!nrow(x)) stop("No terms pass the selected correction/cutoff")
x$plot_adjusted_p <- x[[color_column]]
p <- ggplot(x, aes(x = Count.Overlap.Gene, y = reorder(Gene.Set.Name, Fold.Enrichment),
                   size = Fold.Enrichment, color = plot_adjusted_p)) +
  geom_point() + theme_bw(base_size = 12) +
  labs(title = a$title, x = "Count", y = "DirectGO Term", color = color_column, size = "Fold Enrichment") +
  scale_y_discrete(labels = function(x) vapply(x, function(y) paste(strwrap(y, width = 40), collapse = "\n"), character(1)))
ggsave(a$output, p, device = "pdf", width = 7, height = 7)
write_tsv(x, paste0(a$output, ".plot_input.tsv"))
writeLines(capture.output(dput(a)), paste0(a$output, ".options.R"))
writeLines(capture.output(sessionInfo()), paste0(a$output, ".sessionInfo.txt"))
