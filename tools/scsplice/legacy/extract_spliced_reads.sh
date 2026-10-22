#!/usr/bin/env bash
# Clean ACTIVE samtools.sh step only. Inputs are already strand-separated BAMs.
# CIGAR N means an alignment skips reference sequence (spliced-read selection).
# This does not choose cell barcodes, calculate coverage, or create Figure 6A.
# Usage: bash extract_spliced_reads.sh forward.bam reverse.bam output_prefix THREADS
set -euo pipefail
fail() { printf '%s\n' "$*" >&2; return 1; }
[[ $# == 4 ]] || fail 'Usage: extract_spliced_reads.sh forward.bam reverse.bam output_prefix THREADS'
fw=$1; rv=$2; prefix=$3; threads=$4
[[ $threads =~ ^[1-9][0-9]*$ ]] || fail 'THREADS must be >0'
command -v samtools >/dev/null || fail 'samtools is not available'
[[ -f $fw && -f $rv ]] || fail 'Input BAM missing'
samtools quickcheck "$fw" "$rv"
mkdir -p "$(dirname "$prefix")"
# Use temporary paths so a failed pipeline does not leave an apparent finished BAM.
for strand in fw rv; do
  output="${prefix}.${strand}.bam"
  [[ ! -e $output && ! -e ${output}.bai ]] || fail "Output exists: $output"
done
for strand in fw rv; do
  if [[ $strand == fw ]]; then input=$fw; else input=$rv; fi
  output="${prefix}.${strand}.bam"
  tmp=$(mktemp "${output}.tmp.XXXXXX")
  trap 'rm -f "$tmp" "${tmp}.bai"' ERR INT TERM
  samtools view -@ "$threads" -h "$input" |
    awk 'BEGIN {OFS="\t"} /^@/ {print; next} $6 ~ /N/ {print}' |
    samtools view -@ "$threads" -b -o "$tmp" -
  samtools quickcheck "$tmp"
  samtools index -@ "$threads" "$tmp"
  mv "$tmp" "$output"; mv "${tmp}.bai" "${output}.bai"
  trap - ERR INT TERM
done

