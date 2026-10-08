# Candidate miR-193 target prediction

This workflow predicts miR-193 binding to Heliconius 3′ UTRs and joins gene-level predictions to the comparative differential-expression annotations in Supplementary Table S8.

## Inputs and software

The original inputs and results are in `archive/received/pegasus_miR/`: the mature miRNA FASTA, compressed longest-transcript UTR FASTA, selected transcript IDs, RNAhybrid command, verbose predictions and gene summaries. Python 3.9 or later is sufficient for summarization and annotation. Prediction requires an RNAhybrid environment; version 2.1.2 is reported in the manuscript, but is not established by the supplied job logs.

## Execution

Run from the repository root. Each output directory must be new.

```bash
bash scripts/targets/rnahybrid.sh \
  archive/received/pegasus_miR/HmelFrancGenes.longest.3utr.fa.gz \
  archive/received/pegasus_miR/hme-miR-193.fa \
  results/miR193_predictions.out rnahybrid

python3 scripts/targets/summarize_rnahybrid.py \
  --input archive/received/pegasus_miR/Hmel_miR193_RNAhybrid_longest3utr.out \
  --output-dir results/miR193_summary

python3 scripts/targets/annotate_targets.py \
  --annotations data/paper/Table_S8.tsv \
  --summary results/miR193_summary/gene_hit_summary.tsv \
  --output-dir results/miR193_comparative
```

The launcher uses `micromamba run` to avoid interactive shell activation. It retains `-s 3utr_fly -m 10000 -b 10 -e -25` from the supplied command. Its command record is saved beside the output. The manuscript's −20 kcal/mol description is a distinct setting; no replacement cutoff or additional p-value filter is applied here.

The parser removes the `transcript:` prefix and terminal numeric transcript suffix to recover gene IDs. It exports individual hits, counts per gene, minimum free energy per gene and a predicted-gene list. Empty or incomplete prediction files are rejected.

The comparative annotation script preserves orthology rows and counts each semicolon-separated gene ID once per row. It sums hit counts and takes the lowest energy across mapped genes. The conserved candidate rule is a predicted miR-193 hit plus a DE annotation in at least two of the three studies. Cross-species DE flags are supplied inputs; this step does not repeat those studies' expression analyses.

## Validation and interpretation

Parsing the archived predictions gives 2,079 hits across 1,265 genes, exactly reproducing the supplied summary and predicted-gene set. Comparative annotation reproduces all 733 Table S8 hit counts and Boolean annotations, including the 15 conserved candidate rows. The KrT95D row lacks an energy value in S8; its archived prediction gives −26.4 kcal/mol. Original S8 is unchanged; the comparison output records this difference.

These checks validate result parsing and table joins. RNAhybrid itself has not been rerun here. All three archived job logs contain shell-activation errors and do not establish the successful prediction environment. The original UTR-generation code and primary comparative DE input provenance remain necessary for a complete upstream reconstruction.
