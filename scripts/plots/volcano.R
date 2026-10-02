#!/usr/bin/env Rscript
# Clean volcano plot from saved FindMarkers tables. Vector PDF + exact plotted values.
# Labels TSV (optional): gene_id, label. Optional annotation: gene_id, category.
# Rebuilt plots require comparison with final panel PDFs before claiming exact matches.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages(library(ggplot2))
a <- parse_args(c("input", "p-cutoff", "fc-cutoff", "output")); prepare_output(a$output)
x <- read_tsv(a$input); check_columns(x, c("gene_id", "avg_log2FC", "p_val_adj"))
if (anyDuplicated(x$gene_id)) stop("Duplicate gene IDs")
if (any(!is.finite(x$avg_log2FC)) || any(!is.finite(x$p_val_adj)) || any(x$p_val_adj < 0 | x$p_val_adj > 1)) stop("Invalid DE values")
# Seurat may return p=0 from numeric underflow. Retain originals; floor for display only.
x$plot_neglog10_p <- -log10(pmax(x$p_val_adj, .Machine$double.xmin))
x$category <- ifelse(x$p_val_adj < number(a[["p-cutoff"]], "p-cutoff") & abs(x$avg_log2FC) > number(a[["fc-cutoff"]], "fc-cutoff"), "significant", "other")
if (!is.null(a$annotations)) {
  annotations <- read_tsv(a$annotations); check_columns(annotations, c("gene_id", "category"))
  if (anyDuplicated(annotations$gene_id) || anyNA(annotations$category)) stop("Invalid categories")
  idx <- match(x$gene_id, annotations$gene_id); use <- !is.na(idx)
  x$category[use] <- annotations$category[idx[use]]
}
p <- ggplot(x, aes(avg_log2FC, plot_neglog10_p, color = category)) + geom_point(size = 0.7) +
  theme_classic() + labs(x = "Average log2 fold change (group 1 / group 2)", y = "-log10 adjusted p-value")
if (!is.null(a$labels)) {
  labels <- read_tsv(a$labels); check_columns(labels, c("gene_id", "label"))
  if (anyDuplicated(labels$gene_id)) stop("Duplicate label IDs")
  x$label <- labels$label[match(x$gene_id, labels$gene_id)]
  p <- p + ggrepel::geom_text_repel(data = x[!is.na(x$label), ], aes(label = label), max.overlaps = Inf)
}
ggsave(a$output, p, device = "pdf", width = 7, height = 7)
write_tsv(x, paste0(a$output, ".plot_input.tsv"))
writeLines(capture.output(dput(a)), paste0(a$output, ".options.R"))
writeLines(capture.output(sessionInfo()), paste0(a$output, ".sessionInfo.txt"))
