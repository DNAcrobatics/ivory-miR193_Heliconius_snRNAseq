# scsplice — quick single-cell locus plots

A small standalone extraction of the ivory:mir-193 plotting workflow. Two commands
turn an indexed, barcode-tagged BAM and existing cell annotations into cell-type
read panels and expression ribbons. Extraction scans only the requested locus;
subsequent plotting uses the saved table. No alignment, clustering, annotation
changes, or biological pipeline reruns are needed.

**New to the workflow? Start with [WALKTHROUGH.md](WALKTHROUGH.md)** for a labelled,
step-by-step tutorial covering inputs, barcode matching, Pegasus submission,
reading the plots, zooming and adapting the commands to another locus.

## Install

```bash
micromamba create -f /absolute/path/scsplice/environment.yml
```

Alternatively install `pysam` and `matplotlib` into an existing Python >=3.10
environment. Run the BAM extraction and real-data plotting on Pegasus. The
commands below belong inside a Slurm job; `submit_pegasus.sh` supplies a ready
wrapper. Paths and names below are examples, not claims about existing files.

## Input

One **tab-separated** metadata file with these exact column names:

```tsv
sample	barcode	celltype
iv_20	AAACCCAGTGCGGATA	Scale
iv_20	AAACCCAGTTACCTGA	Socket
iv_33	AAACCCAGTGCGGATA	Scale
```

Sample names must match `--sample`. Barcode matching is exact by default, and is
always scoped to that sample. `--strip-barcode-suffix` removes a terminal `-1`,
`-2`, etc. on both sides only when explicitly requested. Duplicate sample/barcode
pairs (including collisions introduced by stripping) cause an error.

Export existing Seurat annotations if needed, using your actual object and column
names. For the original script's `orig.ident` and `Iv` fields:

```bash
Rscript /absolute/path/scsplice/export_seurat_metadata.R \
  /absolute/path/existing_object.rds /absolute/path/cells.tsv \
  orig.ident Iv rownames last-underscore
```

This exporter requires SeuratObject in the R environment. `last-underscore`
implements the original removal of sample prefixes from Seurat cell names. Use
`exact` when the barcode column/rownames already contain BAM barcodes; use a real
metadata barcode column instead of `rownames` when available. The exporter does
not strip 10x suffixes or modify the object. It has not been executed against the
original Seurat object in this session.

## Two commands

```bash
micromamba run -n scsplice python /absolute/path/scsplice/scsplice.py extract \
  --bam /absolute/path/sample.bam \
  --metadata /absolute/path/cells.tsv --sample iv_20 \
  --region CONTIG:1300000-1700000 \
  --strand forward --strip-barcode-suffix \
  --output /absolute/path/results/iv20

micromamba run -n scsplice python /absolute/path/scsplice/scsplice.py plot \
  --reads /absolute/path/results/iv20.reads.tsv \
  --region CONTIG:1300000-1700000 \
  --annotation /absolute/path/reference.gff3 \
  --output /absolute/path/results/iv20_plot
```

Replace `CONTIG` with the exact contig name in your BAM and annotation. Region
coordinates are **1-based, inclusive**. Annotation is optional and must use the
same reference coordinates. It draws gene lines and exon boxes on separate
forward/reverse rows; it does not infer transcript isoforms.

Add `--spliced-only` to **extract** to retain only alignments with a CIGAR `N`.
Leave it out for the original all-alignment expression ribbons. `--strand`
means SAM alignment orientation (flag 16); it does not infer RNA strand for
paired-end libraries. The original forward plot used a prefiltered forward BAM.

For a ready Pegasus submission (export your own actual paths first):

```bash
export TOOL=/absolute/path/scsplice/scsplice.py
export BAM=/absolute/path/sample.bam
export META=/absolute/path/cells.tsv
export SAMPLE=iv_20
export REGION=CONTIG:1300000-1700000
export PREFIX=/absolute/path/results/iv20
bash /absolute/path/scsplice/submit_pegasus.sh
```

The wrapper uses the `nano` partition, 1 CPU, 8 GB, 30 minutes, the `scsplice`
micromamba environment and alignment orientation `forward`, with suffix stripping
disabled. Override `STRAND`, `STRIP_SUFFIX=1`, `SPLICED_ONLY=1`, `PARTITION`, `MEMORY`,
`TIME`, and `ENV_NAME` as needed. No real job was submitted in this session.

## Replot quickly

Select cell types, zoom, or plot several samples without touching the BAM:

```bash
micromamba run -n scsplice python /absolute/path/scsplice/scsplice.py plot \
  --reads /absolute/path/results/iv20.reads.tsv /absolute/path/results/iv33.reads.tsv \
  --region CONTIG:1400000-1600000 \
  --celltype Scale --celltype Socket \
  --style blocks --binning fixed --bin-width 500 --log-limits 0 3 \
  --output /absolute/path/results/zoom
```

Each read table must represent a different sample. A zoom cannot recover reads
outside the extracted locus. Cell eligibility remains based on the extracted
locus, while row ordering in the plot reflects the visible alignments.

## What the outputs mean

| Output | Meaning |
|---|---|
| `.reads.tsv` | Matched alignment starts, cell types, strands and CIGAR strings |
| `.extract.json` | Extraction options, version and matching/filter counts |
| `.cells.pdf` | One row per cell; solid aligned blocks and dashed CIGAR-N gaps |
| `.ribbons.pdf` | One ribbon per sample × cell type; shared log10 start-count scale |
| `.ribbon.tsv` | Numeric counts underlying the ribbon rectangles |
| `.plot.json` | Plot settings and software version |

`--style starts` uses the original short alignment-start marks instead of aligned
blocks and splice gaps. PDFs use vector elements. Existing output files are
refused; give each iteration a new output prefix.

These are **alignment counts, not UMI counts, per-base coverage, junction usage,
or isoform abundance**. By default cells require at least six matched locus
alignments, reproducing the original `n() > 5`. Change that with `--min-reads`.
No extra duplicate, secondary, supplementary or MAPQ filter is silently added;
`--min-mapq` defaults to zero. Missing CB tags and unmatched reads are counted.

The default `--binning legacy` preserves the source algorithm: sort occupied
alignment-start positions, combine 100 successive occupied positions, sum their
counts and draw a 500 bp rectangle from the first position. Rectangles may overlap
or leave gaps; larger counts are drawn last. It is not uniform genomic binning.
`--binning fixed` explicitly switches to uniform genomic bins, including zeros.
Both methods count alignment starts inside the display interval. Reads overlapping
the interval but starting before it can still appear in the cell panel.

The new tool splits ribbons by cell type and sample; the original loop pooled
the matched cell types within each sample for its summary ribbon. The original
color scale used reversed imola, with log10 limits 2–3. This lightweight version
uses reversed viridis with the same limits; select another Matplotlib palette
with `--palette` and set a suitable range with `--log-limits`. Nonzero values
outside the limits saturate; zero bins keep the light background. This is an
explicit display choice, not a claim of pixel-identical figure reproduction.

## Source and validation

Source repository: https://github.com/Hanliconius/ivory-mir193_paper

Inspected snapshot: `c118c1cec220372846115648194d57f7f08a562f`.
The original script is `archive/drive/Ivory_Snuc/millefy_sc_vis.r`; it contains
barcode export, repeated sample-specific cell panels, read-start ribbons, locus
annotation, `scRibbons.pdf` and an optix-locus extension. The cleaned scripts are
`scripts/locus/plot_barcode_reads.R` and `extract_spliced_reads.sh`.

Unmodified copies are included under `legacy/`; `lib/common.R` supplies the
cleaned R script's relative dependency. Its plotting dependencies include scico,
Rsamtools, GenomicRanges and ggplot2. The old long script is an interactive record,
not a standalone executable. No license has been invented for the source.

The repository's analysis note still marks the final Figure 6A inputs/settings
and result matching as pending. This package is grounded in the accessible
source; it does not establish that every cluster-side file was recovered or
reproduce the final published plot without those inputs.

Validation: synthetic indexed BAM integration test covers repeated barcodes
across samples, explicit suffix matching, strand filtering, splice filtering,
missing CB tags, the six-read threshold, both binning modes, PDF generation and
overwrite refusal. Unit checks verify CIGAR coordinates and historical bins.
No original biological BAM or Seurat object was analyzed locally.

```bash
cd /absolute/path/scsplice
micromamba run -n scsplice python -m unittest -v test_scsplice.py
```
