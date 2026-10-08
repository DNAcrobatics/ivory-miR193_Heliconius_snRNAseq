#!/usr/bin/env python3
"""Summarize verbose RNAhybrid output without adding prediction filters.

Gene identifiers are obtained by removing 'transcript:' and the terminal
numeric transcript suffix. This reproduces the supplied Heliconius summary.
"""
import argparse
import csv
import math
import re
from pathlib import Path


def parse_hits(text):
    blocks = re.split(r'(?m)^target:\s*', text)
    if blocks[0].strip() or len(blocks) == 1:
        raise ValueError('Expected nonempty verbose RNAhybrid target blocks')
    hits = []
    for block in blocks[1:]:
        target = block.splitlines()[0].strip()
        fields = {}
        for key, pattern in [('mfe', r'^mfe:\s*([-\d.eE+]+) kcal/mol\s*$'),
                             ('p_value', r'^p-value:\s*([-\d.eE+]+)\s*$'),
                             ('position', r'^position\s+(\d+)\s*$'),
                             ('mirna', r'^miRNA\s*:\s*(\S+)\s*$')]:
            matches = re.findall(pattern, block, re.M)
            if len(matches) != 1:
                raise ValueError(f'{target}: missing or repeated {key}')
            fields[key] = matches[0]
        if not re.search(r'^target 5', block, re.M) or not re.search(r'^miRNA  3', block, re.M):
            raise ValueError(f'{target}: incomplete alignment')
        mfe, pvalue = float(fields['mfe']), float(fields['p_value'])
        if not math.isfinite(mfe) or not math.isfinite(pvalue) or not 0 <= pvalue <= 1:
            raise ValueError(f'{target}: invalid numeric prediction')
        gene = re.sub(r'\.\d+$', '', target.removeprefix('transcript:'))
        if not re.fullmatch(r'Hmel\.Hmel\w+G\d+', gene):
            raise ValueError(f'Unexpected Heliconius transcript identifier: {target}')
        hits.append(dict(transcript=target, gene_id=gene, mirna=fields['mirna'],
                         mfe=mfe, p_value=pvalue, position=int(fields['position'])))
    return hits


def summarize(hits):
    genes = {}
    for hit in hits:
        entry = genes.setdefault(hit['gene_id'], [0, hit['mfe']])
        entry[0] += 1
        entry[1] = min(entry[1], hit['mfe'])
    return [dict(gene_id=g, n_hits=n, best_MFE=m) for g, (n, m) in sorted(genes.items())]


def write_tsv(path, rows):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    hits = parse_hits(args.input.read_text())
    genes = summarize(hits)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    write_tsv(args.output_dir / 'hits.tsv', hits)
    write_tsv(args.output_dir / 'gene_hit_summary.tsv', genes)
    (args.output_dir / 'predicted_genes.txt').write_text(''.join(x['gene_id'] + '\n' for x in genes))
    print(f'{len(hits)} hits; {len(genes)} genes')


if __name__ == '__main__':
    main()
