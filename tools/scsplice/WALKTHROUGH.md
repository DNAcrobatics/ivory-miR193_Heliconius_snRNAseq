# From a single-cell BAM to labelled locus ribbons

**A practical walkthrough using scsplice**

Suppose we want to inspect transcription across a locus in scale and socket
cells. We already have aligned reads and cell-type assignments. This tutorial
shows how to select a genomic region, connect each read to its cell, and make
two complementary views:

- **Cell panel:** one row per cell, with aligned blocks and dashed splice gaps.
- **Expression ribbon:** one horizontal strip per sample and cell type, coloured
  by the number of alignment starts along the region.

The extraction is done once. We can then change the view without reading the BAM
again. This walkthrough uses Pegasus and submits the analysis through Slurm.

## Before starting: what do we need?

| Input | What it supplies | What to check |
|---|---|---|
| Coordinate-sorted BAM and its index | Aligned reads and cell-barcode tags | Use the BAM for the sample of interest; corrected barcodes are usually in `CB` |
| Cell metadata TSV | The sample and cell type belonging to each barcode | Required columns: `sample`, `barcode`, `celltype` |
| Region | The chromosome/contig and interval to inspect | Must use the same reference assembly as the BAM |
| Optional GFF3/GTF | Gene lines and exon boxes | Contig names and coordinates must match the BAM |

The tool uses existing annotations. It does not recluster cells or realign reads.
The BAM must retain barcode tags; an expression count matrix alone cannot supply
the read positions and splice gaps.

## Step 1 — Put the paths and region in one place

Unpack the toolkit somewhere accessible from Pegasus compute nodes. Edit the
paths below to point to your own files. `/absolute/path/...` values are
placeholders; the sample `iv_20` is an example.

```bash
# Folder containing scsplice.py and submit_pegasus.sh.
export SCSPLICE_DIR=/absolute/path/scsplice

# Existing alignment and cell annotation files.
export BAM=/absolute/path/iv_20.Aligned.sortedByCoord.out.bam
export META=/absolute/path/cells.tsv

# Must match the sample value in cells.tsv exactly.
export SAMPLE=iv_20

# Example from the Hmel2.5 ivory locus: replace for another reference/locus.
export REGION=Hmel215003o:1300000-1700000

# A new output prefix for this extraction. Do not include a file extension.
export PREFIX=/absolute/path/results/ivory_iv20_forward
export TOOL="$SCSPLICE_DIR/scsplice.py"
```

**How to read the region:** `Hmel215003o:1300000-1700000` means bases 1,300,000
through 1,700,000 on contig `Hmel215003o`, including both endpoints. These are
1-based coordinates. The same numbers on another assembly may describe a
different region.

Start from the full aligned BAM when you want a new interval. An older BAM
already cropped to a smaller locus cannot supply reads outside that crop.

## Step 2 — Prepare the barcode-to-cell-type table

The metadata file is a **tab-separated** text file, with a header:

```tsv
sample	barcode	celltype
iv_20	AAACCCAGTGCGGATA	Scale
iv_20	AAACCCAGTTACCTGA	Socket
iv_33	AAACCCAGTGCGGATA	Scale
```

The rows above illustrate the format; they are not real cell assignments. Use
your actual barcodes and labels. The same barcode may occur in different samples.
The tool matches **sample + barcode**, so the `iv_33` row cannot be assigned to an
`iv_20` read.

You can export this table from an existing Seurat RDS object using the included
helper. Skip this section if you already have the TSV. Set the R environment,
object path and metadata column names to your own values:

```bash
export SEURAT_RDS=/absolute/path/annotated_cells.rds
export R_ENV=your_existing_Seurat_environment
export SAMPLE_COLUMN=orig.ident
export CELLTYPE_COLUMN=celltype

mkdir -p "$(dirname "$META")"

sbatch --parsable --partition=nano --cpus-per-task=1 \
  --mem=16G --time=00:30:00 --job-name=cell_metadata \
  --output="${META}.%j.log" --export=ALL <<'SBATCH'
#!/usr/bin/env bash
set -euo pipefail
micromamba run -n "$R_ENV" Rscript "$SCSPLICE_DIR/export_seurat_metadata.R" \
  "$SEURAT_RDS" "$META" \
  "$SAMPLE_COLUMN" "$CELLTYPE_COLUMN" rownames last-underscore
SBATCH
```

Wait for this job to finish before extracting reads. Large Seurat objects may
need more memory or time than the example allocation.

The six exporter arguments, in order, are:

| Argument | Example | Meaning |
|---|---|---|
| Object | `annotated_cells.rds` | An existing Seurat object |
| Output | `cells.tsv` | New metadata file to create |
| Sample column | `orig.ident` | Metadata field identifying the sample |
| Cell-type column | `celltype` | Field containing the labels you want to compare |
| Barcode source | `rownames` | Use cell names; alternatively give a metadata column name |
| Barcode conversion | `last-underscore` | Turn `iv_20_AAACCCAGTGCGGATA` into `AAACCCAGTGCGGATA` |

Use `exact` instead of `last-underscore` if the barcode source already matches
the BAM. For the historical ivory high/low classification, the original script
used the `Iv` column; use that column only when those are the groups you intend
to display. The exporter requires SeuratObject in the selected R environment.

**One barcode detail to get right:** a BAM barcode may be `AAACCCAGTGCGGATA-1`
while the table contains `AAACCCAGTGCGGATA`. Exact matching is the default. Set
`STRIP_SUFFIX=1` in Step 4 only when the terminal `-1`/`-2` suffix should be removed
from both sources. If removal makes two barcodes identical within a sample, the
tool stops rather than merging them.

## Step 3 — Create the plotting environment once

```bash
micromamba create -f "$SCSPLICE_DIR/environment.yml"
```

This creates the Python environment called `scsplice`, containing pysam and
Matplotlib. The optional Seurat exporter uses your separate R environment.

## Step 4 — Submit the first region plot

Choose which reads to retain, then submit:

```bash
# Alignment orientation: forward, reverse, or both.
export STRAND=forward

# Set to 1 only if BAM and metadata differ by a terminal barcode suffix.
export STRIP_SUFFIX=0

# 0 = all matched alignments; 1 = only alignments containing CIGAR N.
export SPLICED_ONLY=0

# The environment created in Step 3.
export ENV_NAME=scsplice

bash "$SCSPLICE_DIR/submit_pegasus.sh"
```

The printed number is the Slurm job ID. The wrapper requests one CPU, 8 GB and
30 minutes on `nano`. It writes the log to `${PREFIX}.JOBID.log` and runs these
two operations in sequence:

1. **Extract:** fetch alignments overlapping the region; apply the chosen
   orientation/splice filters; match sample barcodes to cell types; retain cells
   with at least six matching locus alignments; save a read table.
2. **Plot:** read that table and create a cell panel and ribbon for each cell type.

**Forward and reverse describe alignment orientation**, using SAM flag 16. They
are not an automatic inference of transcript strand for a paired-end library.
Use the orientation appropriate to your library preparation; `both` retains both.

**All alignments versus spliced only:** the historical expression ribbon used
the all-alignment view. Set `SPLICED_ONLY=1` to focus on reads spanning an intron
(`N` in the CIGAR). This also changes which cells pass the six-alignment threshold.
Do not assume an empty splice-only track means the gene is unexpressed.

## Step 5 — Open the outputs and check the match

For the example prefix, the important files are:

| File ending | What to look at |
|---|---|
| `.reads.tsv` | Extracted read records, including barcode, cell type, start and CIGAR |
| `.extract.json` | Counts before/after matching, retained cells, missing barcode tags and options |
| `_plot.cells.pdf` | Each cell is a row; solid blue segments are aligned blocks; dashed grey segments are skipped introns |
| `_plot.ribbons.pdf` | Labelled sample/cell-type strips with a shared colour legend |
| `_plot.ribbon.tsv` | The numeric counts used to draw each strip |
| `_plot.plot.json` | Settings used to create the figure |

Both plots label the genomic axis. The ribbon legend reports **log10 alignment
starts**. A value of 2 corresponds to 100 starts; a value of 3 corresponds to
1,000 starts in the displayed ribbon rectangle.

Read `.extract.json` before interpreting an unexpectedly sparse panel. Many
unmatched reads can indicate a metadata/sample mismatch, while missing `CB` tags
can indicate the wrong BAM or barcode tag. The saved counts are for this region
and the requested filters, not the whole library.

The cell count on each panel includes cells with visible alignments after plotting
filters. Cells with no retained locus reads do not appear as blank rows.

## Step 6 — Zoom in or compare selected cell types

We now have the reads we need. Submit a plotting-only job to make a closer view:

```bash
# A smaller interval inside the extracted region.
export VIEW_REGION=Hmel215003o:1400000-1600000
export VIEW_PREFIX="${PREFIX}_zoom"

sbatch --parsable --partition=nano --cpus-per-task=1 \
  --mem=8G --time=00:30:00 --job-name=splice_zoom \
  --output="${VIEW_PREFIX}.%j.log" --export=ALL <<'SBATCH'
#!/usr/bin/env bash
set -euo pipefail
micromamba run -n "$ENV_NAME" python "$TOOL" plot \
  --reads "${PREFIX}.reads.tsv" \
  --region "$VIEW_REGION" \
  --celltype Scale --celltype Socket \
  --style blocks \
  --output "$VIEW_PREFIX"
SBATCH
```

Replace `Scale` and `Socket` with labels present in your table; spelling and case
must match. Remove both `--celltype` options to show all cell types. Repeat the
option to include more groups. A new output prefix keeps the first plot intact.

The underlying plotting command also accepts these useful options. Add them
before the final `--output` line in the job above:

| Add to `plot` | Effect |
|---|---|
| `--style starts` | Show short alignment-start marks, as in the historical panel |
| `--annotation /absolute/path/reference.gff3` | Add gene/exon context on forward/reverse rows |
| `--log-limits 0 3` | Use a shared colour scale spanning 1–1,000 starts |
| `--binning fixed --bin-width 500` | Use uniform 500 bp genomic bins |
| `--reads /path/sample1.reads.tsv /path/sample2.reads.tsv` | Stack samples; each table must represent a different sample |

For multiple samples, run extraction once per sample with its matching BAM,
sample ID and a distinct prefix, preferably using the same region and filters.
Then supply the resulting tables together to `plot`.

Zooming uses the saved extraction: it cannot recover reads beyond that original
interval. The six-read eligibility threshold is still based on the extraction
region. Extract again with a new prefix if you want eligibility defined by a
different interval, orientation or splice filter.

## Step 7 — Adapt the workflow to another gene

Keep the same BAM and metadata if you are looking at another locus in the same
sample. Change the region and output prefix, then submit again:

```bash
# Replace these with the actual coordinates and a new output name.
export REGION=YOUR_CONTIG:100000-200000
export PREFIX=/absolute/path/results/my_next_gene

bash "$SCSPLICE_DIR/submit_pegasus.sh"
```

The wrapper uses the current `STRAND`, `STRIP_SUFFIX` and `SPLICED_ONLY` settings.
Set them explicitly when the next locus requires different choices. If changing
samples, update `BAM` and `SAMPLE` too.

## How to interpret the ribbons

The ribbon summarises **alignment starts**, rather than per-base coverage or
UMI expression. Its intensity depends on the number of included cells, sequencing
depth and filtering. The shared colour scale makes the values readable across
tracks, but does not normalise those differences. Use it to inspect where reads
accumulate; it is not a differential-expression or isoform-usage test.

Two binning modes are available:

| Mode | Calculation | When to use it |
|---|---|---|
| `legacy` — default | Combine 100 successive occupied start positions; draw a 500 bp rectangle from the first position | Retain the historical ivory ribbon calculation |
| `fixed` | Count starts in consecutive fixed-width genomic bins | Make the relationship between each rectangle and its genomic interval explicit |

Legacy rectangles can overlap or leave gaps, and their 500 bp drawing width does
not delimit the positions that contributed to the count. Choose one method
consistently when comparing tracks. Both methods count starts inside the display
interval; an overlapping alignment that starts outside it can appear in the cell
panel without contributing a ribbon start.

The default colour limits are log10 2–3. Nonzero values outside the chosen range
saturate at the end colours; zero bins keep the pale background. For sparse
regions, try `--log-limits 0 3` and inspect the numeric `.ribbon.tsv` alongside
the plot. The new tool uses reversed viridis; the original used reversed imola.

## If something does not look right

| Symptom | First thing to check |
|---|---|
| Region does not match the BAM | Assembly, exact contig name and contig length |
| No metadata for this sample | `SAMPLE` must match the `sample` column |
| No reads survive matching | Barcode format/suffixes, orientation, splice selection and six-read minimum |
| Duplicate barcode error | Repeated rows or a suffix-removal collision within a sample |
| Missing barcode tags | Whether the BAM retains corrected `CB` tags; `extract --barcode-tag TAG` can select another tag explicitly |
| Requested cell type is absent | Exact label and whether that group has retained reads in the region |
| Ribbons show little colour variation | Numeric counts versus the selected log colour limits |
| Output already exists | Choose a new prefix for the extraction or plot |
| Job exceeds memory/time | Inspect the Slurm log; override `MEMORY`, `TIME` and `PARTITION` for the submission wrapper |

For the wrapper, resource overrides are environment variables, for example
`export MEMORY=16G`. For manual extraction or non-default thresholds, inspect
`python scsplice.py extract --help`; place the resulting command in a Slurm job
like the plotting example above. `--min-reads` defaults to 6 and `--min-mapq` to 0.
No duplicate, secondary or supplementary alignment removal is added implicitly.

## What has been checked?

The tool has passed synthetic tests for sample-specific barcode matching,
suffix handling, orientation/splice selection, CIGAR coordinates, both binning
modes and PDF output. The R exporter has not been run against the original Seurat
object. Reproducing the final ivory figure still requires the original biological
inputs and settings; this walkthrough does not claim that reproduction.

See [README.md](README.md) for the full command reference and source provenance,
and [VALIDATION.md](VALIDATION.md) for the checks already performed.
