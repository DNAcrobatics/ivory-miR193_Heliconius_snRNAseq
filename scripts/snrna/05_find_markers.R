#!/usr/bin/env Rscript
# Clean FindMarkers calls with two explicit identities and a frozen membership table.
# This implements nucleus-level Wilcoxon DE, not replicate-level pseudobulk DE.
# Membership TSV: barcode and group; group values must equal --group1 / --group2.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages(library(Seurat))
a <- parse_args(c("input", "membership", "assay", "group1", "group2", "min-pct", "logfc-threshold", "output"))
prepare_output(a$output); obj <- read_object(a$input)
m <- read_tsv(a$membership); check_columns(m, c("barcode", "group"))
if (!nrow(m) || anyNA(m[c("barcode", "group")]) || anyDuplicated(m$barcode) || any(!m$barcode %in% colnames(obj))) stop("Invalid cell memberships")
if (a$group1 == a$group2 || !setequal(unique(m$group), c(a$group1, a$group2))) stop("Membership must contain exactly the two requested groups")
if (any(table(m$group) < 3L)) stop("Fewer than three nuclei in one group")
obj <- subset(obj, cells = m$barcode)
if (!a$assay %in% Assays(obj)) stop("Assay missing")
DefaultAssay(obj) <- a$assay
# Layer joining is deliberately left to verified object preparation: do not merge layers silently.
obj$comparison_group <- m$group[match(colnames(obj), m$barcode)]
Idents(obj) <- "comparison_group"
markers <- FindMarkers(obj, ident.1 = a$group1, ident.2 = a$group2,
  assay = a$assay, test.use = "wilcox", recorrect_umi = FALSE,
  min.pct = number(a[["min-pct"]], "min-pct"),
  logfc.threshold = number(a[["logfc-threshold"]], "logfc-threshold"))
# Positive fold changes mean group1 > group2. Keep gene identifiers as an explicit column.
write_tsv(cbind(gene_id = rownames(markers), markers), a$output)
write_tsv(cbind(barcode = colnames(obj), obj[[]]), paste0(a$output, ".membership.tsv"))
writeLines(capture.output(dput(a)), paste0(a$output, ".options.R"))
writeLines(capture.output(sessionInfo()), paste0(a$output, ".sessionInfo.txt"))
