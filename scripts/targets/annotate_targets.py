#!/usr/bin/env python3
"""Join RNAhybrid gene summaries to the comparative DE annotations in Table S8.

Cross-species DE flags are inputs, not recalculated tests. Multiple gene IDs
within one orthology row are counted once each. Conserved candidates have a
predicted hit and DE annotations in at least two of the three studies.
"""
import argparse
import csv
from pathlib import Path

FLAGS = ['Differentialy expressed in this study?', 'Differentially expressed in Tian et al?',
         'Differentially expressed in Fandino et al?']
DERIVED = ['n_RNAhybrid_hits', 'best_RNAhybrid_MFE', 'miR193_RNAhybrid_hit',
           'shared_DE_count', 'shared_DE', 'shared_all_three', 'candidate_conserved_miR193_target']


def boolean(value):
    if value.lower() not in ('true', 'false'):
        raise ValueError(f'Expected True/False DE annotation, received {value!r}')
    return value.lower() == 'true'


def annotate(row, summary):
    row = dict(row)
    genes = set(x.strip() for x in row['gene'].split(';') if x.strip())
    matches = [summary[g] for g in genes if g in summary]
    count = sum(int(x['n_hits']) for x in matches)
    shared = sum(boolean(row[x]) for x in FLAGS)
    row.update(n_RNAhybrid_hits=count,
               best_RNAhybrid_MFE=min(float(x['best_MFE']) for x in matches) if matches else '',
               miR193_RNAhybrid_hit=count > 0, shared_DE_count=shared,
               shared_DE=shared >= 2, shared_all_three=shared == 3,
               candidate_conserved_miR193_target=shared >= 2 and count > 0)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotations', type=Path, required=True)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    with args.summary.open() as stream:
        summary = {r['gene_id']: r for r in csv.DictReader(stream, delimiter='\t')}
    with args.annotations.open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        required = {'gene', *FLAGS}
        if not required.issubset(reader.fieldnames):
            raise ValueError(f'Missing columns: {required - set(reader.fieldnames)}')
        originals = list(reader)
        fields = list(reader.fieldnames) + [x for x in DERIVED if x not in reader.fieldnames]
    rows = [annotate(r, summary) for r in originals]
    differences = []
    for number, (old, new) in enumerate(zip(originals, rows), start=1):
        for field in DERIVED:
            if field not in old:
                continue
            a, b = old[field], str(new[field])
            equal = a.lower() == b.lower() or (a.lower() in ('', 'na', 'nan') and b == '')
            if not equal and a and b:
                try:
                    equal = float(a) == float(b)
                except ValueError:
                    pass
            if not equal:
                differences.append(dict(row=number, gene=old['gene'], field=field, supplied=a, derived=b))
    args.output_dir.mkdir(parents=True, exist_ok=False)
    for filename, selected in [('annotated_targets.tsv', rows), ('conserved_candidates.tsv',
                               [r for r in rows if r['candidate_conserved_miR193_target']])]:
        with (args.output_dir / filename).open('x', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
            writer.writeheader()
            writer.writerows(selected)
    with (args.output_dir / 'annotation_comparison.tsv').open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['row', 'gene', 'field', 'supplied', 'derived'], delimiter='\t')
        writer.writeheader()
        writer.writerows(differences)
    print(f'{len(rows)} rows; {sum(r["candidate_conserved_miR193_target"] for r in rows)} candidates; {len(differences)} annotation differences')


if __name__ == '__main__':
    main()
