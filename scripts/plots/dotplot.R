#!/usr/bin/env Rscript
# Clean grouped Seurat DotPlot; gene identities and display labels stay separate.
# Membership: barcode, group. Gene table: gene_id, label in the intended plot order.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
a <- parse_args(c("input", "membership", "genes", "assay", "output")); prepare_output(a$output)
obj <- read_object(a$input); m <- read_tsv(a$membership); g <- read_tsv(a$genes)
check_columns(m, c("barcode", "group")); check_columns(g, c("gene_id", "label"))
if (!nrow(m) || anyNA(m) || anyDuplicated(m$barcode) || any(!m$barcode %in% colnames(obj))) stop("Invalid membership")
if (!nrow(g) || anyNA(g) || anyDuplicated(g$gene_id)) stop("Invalid gene table")
if (!a$assay %in% Assays(obj)) stop("Assay missing")
DefaultAssay(obj) <- a$assay
if (any(!g$gene_id %in% rownames(obj))) stop("Requested genes absent from assay")
obj <- subset(obj, cells = m$barcode)
obj$comparison_group <- m$group[match(colnames(obj), m$barcode)]
Idents(obj) <- factor(obj$comparison_group, levels = unique(m$group))
p <- DotPlot(obj, features = g$gene_id, assay = a$assay) +
  scale_x_discrete(labels = setNames(g$label, g$gene_id)) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(a$output, p, device = "pdf", width = 7, height = 7)
write_tsv(p$data, paste0(a$output, ".plot_input.tsv"))
writeLines(capture.output(dput(a)), paste0(a$output, ".options.R"))
writeLines(capture.output(sessionInfo()), paste0(a$output, ".sessionInfo.txt"))
