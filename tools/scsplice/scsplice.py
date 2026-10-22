#!/usr/bin/env python3
"""Cell-type read panels and read-start ribbons derived from the ivory workflow."""
import argparse
import collections
import csv
import json
import math
from pathlib import Path
import re
import sys

VERSION = '0.1.0'


def table(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_table(path, rows, columns):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def new_paths(prefix, suffixes):
    paths = [Path(str(prefix) + s) for s in suffixes]
    for p in paths:
        if p.exists():
            raise ValueError(f'Output exists: {p}; choose a new prefix')
    paths[0].parent.mkdir(parents=True, exist_ok=True)
    return paths


def region(value):
    match = re.fullmatch(r'(.+):(\d+)-(\d+)', value.replace(',', ''))
    if not match or not 1 <= int(match[2]) <= int(match[3]):
        raise argparse.ArgumentTypeError('Use contig:start-end (1-based, inclusive)')
    return match[1], int(match[2]), int(match[3])


def barcode(value, strip):
    return re.sub(r'-\d+$', '', value) if strip else value


def extract(a):
    import pysam
    contig, start, end = a.region
    meta = table(a.metadata)
    required = {'sample', 'barcode', 'celltype'}
    if not meta or not required <= set(meta[0]):
        raise ValueError('Metadata requires tab-separated sample, barcode, celltype columns')
    selected = {}
    for r in meta:
        if r['sample'] != a.sample:
            continue
        bc = barcode(r['barcode'], a.strip_barcode_suffix)
        if not bc or bc == '-' or not r['celltype']:
            raise ValueError('Empty barcode or cell type in selected sample')
        if bc in selected:
            raise ValueError(f'Duplicate barcode in sample after matching rule: {bc}')
        selected[bc] = r['celltype']
    if not selected:
        raise ValueError('No metadata rows match --sample')
    reads = []
    audit = collections.Counter()
    with pysam.AlignmentFile(a.bam, 'rb') as bam:
        if not bam.has_index():
            raise ValueError('BAM needs an index')
        if contig not in bam.references or end > bam.get_reference_length(contig):
            raise ValueError('Region does not match BAM reference names/lengths')
        for r in bam.fetch(contig, start - 1, end):
            audit['region_alignments'] += 1
            if r.is_unmapped:
                continue
            if a.strand != 'both' and r.is_reverse != (a.strand == 'reverse'):
                continue
            if r.mapping_quality < a.min_mapq:
                continue
            if a.spliced_only and not any(op == 3 for op, n in (r.cigartuples or [])):
                continue
            audit['filtered_alignments'] += 1
            if not r.has_tag(a.barcode_tag):
                audit['missing_barcode_tag'] += 1
                continue
            bc = barcode(str(r.get_tag(a.barcode_tag)), a.strip_barcode_suffix)
            if bc not in selected:
                audit['unmatched_alignments'] += 1
                continue
            reads.append(dict(sample=a.sample, barcode=bc, celltype=selected[bc],
                              contig=contig, pos=r.reference_start + 1,
                              cigar=r.cigarstring, strand='-' if r.is_reverse else '+'))
    counts = collections.Counter(r['barcode'] for r in reads)
    reads = [r for r in reads if counts[r['barcode']] >= a.min_reads]
    if not reads:
        raise ValueError('No reads survive metadata matching and --min-reads; inspect sample/barcodes/filters')
    paths = new_paths(a.output, ['.reads.tsv', '.extract.json'])
    write_table(paths[0], reads, ['sample', 'barcode', 'celltype', 'contig', 'pos', 'cigar', 'strand'])
    audit.update(retained_alignments=len(reads), retained_cells=len(set(r['barcode'] for r in reads)),
                 metadata_cells=len(selected))
    run = {'version': VERSION, 'options': vars(a), 'counts': dict(audit),
           'pysam_version': pysam.__version__,
           'coordinates': '1-based inclusive input and alignment start; BAM fetch includes overlapping reads',
           'metric': 'alignments, not UMIs; no implicit duplicate/secondary/supplementary removal'}
    paths[1].write_text(json.dumps(run, indent=2, default=str) + '\n')
    print(f'{len(reads)} alignments; {audit["retained_cells"]} cells -> {paths[0]}')


def ribbon_bins(reads, method, bin_size, width, start, end):
    counts = collections.Counter(int(r['pos']) for r in reads if start <= int(r['pos']) <= end)
    if method == 'legacy':
        positions = sorted(counts)
        bins = [(positions[i], positions[i] + width, sum(counts[p] for p in positions[i:i+bin_size]))
                for i in range(0, len(positions), bin_size)]
    else:
        bins = []
        for left in range(start, end + 1, width):
            bins.append((left, min(left + width, end + 1), 0))
        totals = collections.Counter()
        for p, n in counts.items():
            totals[(p - start) // width] += n
        bins = [(l, r, totals[i]) for i, (l, r, _) in enumerate(bins)]
    return sorted(bins, key=lambda b: b[2])


def cigar_segments(pos, cigar):
    """Return aligned blocks and skipped introns, as 0-based half-open intervals."""
    x = int(pos) - 1
    blocks, introns = [], []
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        n = int(n)
        if op in 'M=X':
            blocks.append((x, x + n))
            x += n
        elif op == 'N':
            introns.append((x, x + n))
            x += n
        elif op == 'D':
            x += n
    return blocks, introns


def plot(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.patches import Rectangle
    from matplotlib.colors import Normalize
    contig, start, end = a.region
    rows = []
    seen_samples = set()
    for path in a.reads:
        records = table(path)
        samples = {r.get('sample') for r in records}
        if seen_samples & samples:
            raise ValueError('Repeated sample across read tables; supply each sample only once')
        seen_samples.update(samples)
        rows.extend(records)
    required = {'sample', 'barcode', 'celltype', 'contig', 'pos', 'cigar'}
    if not rows or not required <= set(rows[0]):
        raise ValueError('Use a .reads.tsv produced by extract')
    # The prepared table defines the locus used for cell eligibility; plotting can zoom.
    rows = [r for r in rows if r['contig'] == contig and
            int(r['pos']) <= end and any(right >= start for left, right in
                                         cigar_segments(r['pos'], r['cigar'])[0])]
    if a.celltype:
        rows = [r for r in rows if r['celltype'] in a.celltype]
    if not rows:
        raise ValueError('No alignments for selected region/cell types')
    groups = collections.OrderedDict()
    for r in rows:
        key = (r['sample'], r['celltype'])
        groups.setdefault(key, []).append(r)
    paths = new_paths(a.output, ['.ribbons.pdf', '.cells.pdf', '.ribbon.tsv', '.plot.json'])
    cmap = plt.get_cmap(a.palette)
    norm = Normalize(*a.log_limits, clip=True)
    fig, axes = plt.subplots(len(groups), 1, figsize=(a.width, max(2, len(groups)*0.65)),
                             sharex=True, squeeze=False)
    bin_rows = []
    for ax, ((sample, celltype), rs) in zip(axes[:, 0], groups.items()):
        bins = ribbon_bins(rs, a.binning, a.occupied_positions, a.bin_width, start, end)
        ax.set_facecolor('#f5f5f5')
        patches, colors = [], []
        for left, right, count in bins:
            bin_rows.append(dict(sample=sample, celltype=celltype, start=left, end_exclusive=right, count=count))
            if count:
                patches.append(Rectangle((left, 0), right-left, 1))
                colors.append(cmap(norm(math.log10(count))))
        ax.add_collection(PatchCollection(patches, facecolor=colors, edgecolor='none'))
        ax.set(xlim=(start, end+1), ylim=(0, 1), yticks=[], ylabel=f'{sample}\n{celltype}')
        ax.spines[['top', 'right', 'left']].set_visible(False)
    axes[-1, 0].set_xlabel(f'{contig} position (bp) — alignment starts')
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=list(axes[:, 0]),
                 label='log10 alignment starts', fraction=0.025, pad=0.03)
    fig.savefig(paths[0], bbox_inches='tight')
    plt.close(fig)
    heights = [max(1, min(6, len({r['barcode'] for r in rs}) / 45)) for rs in groups.values()]
    fig, axes = plt.subplots(len(groups) + bool(a.annotation), 1, figsize=(a.width, sum(heights)+1.5),
                             sharex=True, squeeze=False,
                             gridspec_kw={'height_ratios': heights + ([0.7] if a.annotation else [])})
    for ax, ((sample, celltype), rs) in zip(axes[:, 0], groups.items()):
        counts = collections.Counter(r['barcode'] for r in rs)
        cells = sorted(counts, key=lambda bc: (counts[bc], bc))
        levels = {bc: i for i, bc in enumerate(cells)}
        solid, skipped, starts, ys = [], [], [], []
        for r in rs:
            y = levels[r['barcode']]
            if a.style == 'starts':
                starts.append(int(r['pos']))
                ys.append(y)
            else:
                blocks, introns = cigar_segments(r['pos'], r['cigar'])
                solid.extend([[(l+1, y), (h, y)] for l, h in blocks])
                skipped.extend([[(l+1, y), (h, y)] for l, h in introns])
        if a.style == 'blocks':
            ax.add_collection(LineCollection(skipped, colors='#999999', linewidths=0.35, linestyles='dashed'))
            ax.add_collection(LineCollection(solid, colors='#2354a4', linewidths=0.8))
        else:
            ax.scatter(starts, ys, marker='_', s=2, color='#2354a4', linewidths=0.5)
        ax.set(xlim=(start, end+1), ylim=(-1, len(cells)), yticks=[],
               ylabel=f'{sample}\n{celltype}\n{len(cells)} cells')
        ax.spines[['top', 'right', 'left']].set_visible(False)
    if a.annotation:
        ax = axes[-1, 0]
        with open(a.annotation) as f:
            for line in f:
                if line.startswith('##FASTA'):
                    break
                if line.startswith('#') or not line.strip():
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) != 9 or cols[0] != contig:
                    continue
                left, right = int(cols[3]), int(cols[4])
                if right < start or left > end:
                    continue
                y = 0.75 if cols[6] == '+' else 0.25
                if cols[2] == 'gene':
                    ax.plot([left, right], [y, y], color='black', lw=0.6)
                elif cols[2] == 'exon':
                    ax.add_patch(Rectangle((left, y-0.1), right-left+1, 0.2, color='black'))
        ax.set(ylim=(0, 1), yticks=[0.25, 0.75], yticklabels=['−', '+'], ylabel='Annotation')
        ax.spines[['top', 'right', 'left']].set_visible(False)
    axes[-1, 0].set_xlabel(f'{contig} position (bp)')
    fig.tight_layout()
    fig.savefig(paths[1], bbox_inches='tight')
    plt.close(fig)
    write_table(paths[2], bin_rows, ['sample', 'celltype', 'start', 'end_exclusive', 'count'])
    paths[3].write_text(json.dumps({'version': VERSION, 'options': vars(a),
                                  'matplotlib_version': matplotlib.__version__,
                                  'metric': 'alignment starts; ribbons do not count junctions or UMIs'},
                                 indent=2, default=str) + '\n')
    print(f'{len(groups)} tracks -> {paths[0]} and {paths[1]}')


def positive(s):
    n = int(s)
    if n < 1:
        raise argparse.ArgumentTypeError('Must be >=1')
    return n


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='version', version=VERSION)
    sub = p.add_subparsers(dest='command', required=True)
    e = sub.add_parser('extract', help='Read one indexed BAM and match cells within one sample')
    for key in ['bam', 'metadata', 'sample', 'output']:
        e.add_argument('--'+key, required=True)
    e.add_argument('--region', type=region, required=True)
    e.add_argument('--strand', choices=['forward', 'reverse', 'both'], default='both',
                   help='Alignment orientation (SAM flag 16), not inferred transcript strand')
    e.add_argument('--spliced-only', action='store_true', help='Retain CIGAR N alignments only')
    e.add_argument('--barcode-tag', default='CB')
    e.add_argument('--strip-barcode-suffix', action='store_true', help='Explicitly remove terminal -digits on both sides')
    e.add_argument('--min-reads', type=positive, default=6, help='Minimum matched locus alignments per cell (default 6 = historical >5)')
    e.add_argument('--min-mapq', type=int, default=0)
    q = sub.add_parser('plot', help='Replot saved reads without rescanning BAMs')
    q.add_argument('--reads', required=True, nargs='+', help='One extracted table per sample')
    q.add_argument('--region', type=region, required=True)
    q.add_argument('--output', required=True)
    q.add_argument('--celltype', action='append', help='Repeat to select multiple cell types')
    q.add_argument('--style', choices=['starts', 'blocks'], default='blocks')
    q.add_argument('--binning', choices=['legacy', 'fixed'], default='legacy')
    q.add_argument('--occupied-positions', type=positive, default=100)
    q.add_argument('--bin-width', type=positive, default=500)
    q.add_argument('--log-limits', nargs=2, type=float, default=[2, 3])
    q.add_argument('--palette', default='viridis_r', help='Matplotlib palette; original used reversed imola')
    q.add_argument('--annotation', help='GFF3/GTF gene and exon features on the same reference')
    q.add_argument('--width', type=float, default=12)
    a = p.parse_args()
    try:
        if a.command == 'plot' and (a.width <= 0 or not a.log_limits[0] < a.log_limits[1]):
            raise ValueError('Require width >0 and increasing --log-limits')
        if a.command == 'extract' and not 0 <= a.min_mapq <= 255:
            raise ValueError('Require 0 <= --min-mapq <= 255')
        (extract if a.command == 'extract' else plot)(a)
    except (ValueError, OSError, ImportError) as exc:
        p.exit(2, f'error: {exc}\n')


if __name__ == '__main__':
    main()
