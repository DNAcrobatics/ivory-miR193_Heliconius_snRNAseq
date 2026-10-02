#!/usr/bin/env Rscript
# Clean the S8 AggregateExpression -> CPM -> ranking/histogram branch.
# Supply the intended all-cell or WT-only object explicitly; no genotype inference.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
a <- parse_args(c("input", "group-by", "group", "output")); prepare_output(a$output)
obj <- read_object(a$input); check_columns(obj[[]], a[["group-by"]])
DefaultAssay(obj) <- "RNA"
counts <- AggregateExpression(obj, assays = "RNA", group.by = a[["group-by"]], slot = "counts")$RNA
cpm <- edgeR::cpm(counts)
if (!a$group %in% colnames(cpm)) stop("Group missing in aggregated matrix; inspect group naming")
values <- sort(cpm[, a$group], decreasing = TRUE); values <- values[values > 1]
tab <- data.frame(gene_id = names(values), cpm = as.numeric(values), rank = seq_along(values))
write_tsv(tab, a$output)
if (!"Ivory" %in% tab$gene_id) stop("Ivory absent from >1 CPM ranking; table saved")
p <- ggplot(tab, aes(cpm)) + geom_histogram(bins = 30) + scale_x_log10() +
  geom_vline(xintercept = tab$cpm[tab$gene_id == "Ivory"]) +
  labs(x = paste("Aggregate expression in", a$group, "(CPM)"), y = "Genes") + theme_classic()
ggsave(paste0(a$output, ".pdf"), p, width = 5, height = 2, device = "pdf")
writeLines(capture.output(dput(a)), paste0(a$output, ".options.R"))
writeLines(capture.output(sessionInfo()), paste0(a$output, ".sessionInfo.txt"))
