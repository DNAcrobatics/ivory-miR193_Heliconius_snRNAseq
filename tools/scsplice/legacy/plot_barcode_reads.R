#!/usr/bin/env Rscript
# Reconstructed from the loop in archive/drive/Ivory_Snuc/millefy_sc_vis.r.
# Plots alignment STARTS, not basewise coverage or UMI counts. Supply a strand-
# filtered, coordinate-sorted/indexed BAM containing CB tags. No additional
# mapping-quality, duplicate or splicing filter is introduced here.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script), "..", "lib", "common.R"))
suppressPackageStartupMessages({library(Rsamtools); library(ggplot2); library(scico)})
a <- parse_args(c("bam", "metadata", "sample", "contig", "start", "end", "output-prefix"))
start <- integer_option(a$start, "start"); end <- integer_option(a$end, "end")
if (end < start) stop("end must be >= start")
prefix <- a[["output-prefix"]]
for (suffix in c(".reads.pdf", ".ribbon.pdf", ".read_starts.tsv", ".ribbon.tsv", ".audit.tsv", ".options.R", ".sessionInfo.txt")) prepare_output(paste0(prefix, suffix))
# bcdf.all.txt is whitespace separated and stores bare 16-base barcodes.
# Restrict to the sample BEFORE joining, since barcodes recur across libraries.
meta <- read.table(a$metadata, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
check_columns(meta, c("barcode", "indiv", "celltype"))
meta <- meta[meta$indiv == a$sample, , drop = FALSE]
if (!nrow(meta)) stop("No metadata for sample")
meta$barcode <- sub("-[0-9]+$", "", meta$barcode)
if (anyNA(meta$barcode) || anyDuplicated(meta$barcode)) stop("Metadata must have unique barcodes within sample")
region <- GenomicRanges::GRanges(a$contig, IRanges::IRanges(start, end))
b <- scanBam(a$bam, param = ScanBamParam(which = region, what = "pos", tag = "CB"))
if (any(vapply(b, function(x) length(x$tag$CB) != length(x$pos), logical(1)))) stop("BAM is missing CB tags")
reads <- do.call(rbind, lapply(b, function(x) data.frame(pos = x$pos, barcode = x$tag$CB, stringsAsFactors = FALSE)))
raw_n <- nrow(reads)
reads <- reads[!is.na(reads$barcode) & !reads$barcode %in% c("-", ""), , drop = FALSE]
reads$barcode <- sub("-[0-9]+$", "", reads$barcode)
# The source retains barcodes with >5 alignments in the input locus BAM.
# Use the supplied interval to define that input locus explicitly.
counts <- table(reads$barcode)
reads <- reads[reads$barcode %in% names(counts)[counts > 5], , drop = FALSE]
eligible_n <- nrow(reads)
reads <- merge(reads, meta, by = "barcode", sort = FALSE)
if (!nrow(reads)) stop("No eligible reads matched this sample's metadata")
counts <- table(reads$barcode)
reads$barcode <- factor(reads$barcode, levels = names(sort(counts, decreasing = FALSE)))
p <- ggplot(reads, aes(pos, barcode)) + geom_point(size = 0.3, color = "blue", shape = 95) +
  facet_grid(celltype ~ ., scales = "free", space = "free_y") + coord_cartesian(xlim = c(start, end)) +
  theme_classic() + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
  axis.ticks.y = element_blank(), strip.background = element_blank()) +
  labs(x = paste0(a$contig, " alignment start"), title = a$sample)
ggsave(paste0(prefix, ".reads.pdf"), p, device = "pdf", width = 20, height = 3)
# Historical ribbons combine 100 successive OCCUPIED start positions, then
# draw a 500 bp rectangle. These are not uniform 100 bp genomic bins.
positions <- as.data.frame(table(reads$pos), stringsAsFactors = FALSE)
positions$pos <- as.numeric(positions$Var1)
positions <- positions[order(positions$pos), ]
groups <- ceiling(seq_len(nrow(positions)) / 100)
bins <- do.call(rbind, lapply(split(positions, groups), function(x)
  data.frame(pos = x$pos[1], count = sum(x$Freq))))
bins <- bins[order(bins$count), ]
r <- ggplot(bins, aes(xmin = pos, xmax = pos + 500, ymin = 0.975, ymax = 1.025, fill = log10(count))) +
  annotate("rect", xmin = start, xmax = end, ymin = 0.975, ymax = 1.025, fill = "#F5F5F5") +
  geom_rect() + scale_fill_gradientn(colors = scico::scico(256, palette = "imola", direction = -1), limits = c(2, 3)) +
  coord_cartesian(xlim = c(start, end)) + theme_void() + theme(legend.position = "none")
ggsave(paste0(prefix, ".ribbon.pdf"), r, device = "pdf", width = 20, height = 1)
write_tsv(reads, paste0(prefix, ".read_starts.tsv"))
write_tsv(bins, paste0(prefix, ".ribbon.tsv"))
write_tsv(data.frame(sample = a$sample, input_alignments = raw_n, eligible_alignments = eligible_n,
                    matched_alignments = nrow(reads), matched_barcodes = length(unique(reads$barcode))), paste0(prefix, ".audit.tsv"))
writeLines(capture.output(dput(a)), paste0(prefix, ".options.R"))
writeLines(capture.output(sessionInfo()), paste0(prefix, ".sessionInfo.txt"))

