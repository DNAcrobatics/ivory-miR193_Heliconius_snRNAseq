#!/usr/bin/env bash
# Predict miR-193 binding with the parameters in the supplied Pegasus command.
# Run within the recorded RNAhybrid environment. A compressed target is accepted.
# This launcher fixes shell activation; it does not change prediction parameters.
set -euo pipefail
if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo 'Usage: rnahybrid.sh targets.fa[.gz] mirna.fa output.out [micromamba-env]' >&2
  exit 2
fi
target=$1
query=$2
output=$3
environment=${4:-rnahybrid}
[[ -s "$target" && -s "$query" ]] || { echo 'Missing or empty FASTA input' >&2; exit 2; }
[[ ! -e "$output" && ! -e "$output.command.txt" ]] || { echo 'Output already exists' >&2; exit 2; }
command -v micromamba >/dev/null
mkdir -p "$(dirname "$output")"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
if [[ "$target" == *.gz ]]; then
  gzip -dc -- "$target" > "$tmp/targets.fa"
  target="$tmp/targets.fa"
fi
command=(micromamba run -n "$environment" RNAhybrid -t "$target" -q "$query" -s 3utr_fly -m 10000 -b 10 -e -25)
printf '%q ' "${command[@]}" > "$output.command.txt"
printf '\n' >> "$output.command.txt"
# The destination is published only after the process exits successfully.
"${command[@]}" > "$tmp/predictions.out"
[[ -s "$tmp/predictions.out" ]] || { echo 'RNAhybrid produced no output' >&2; exit 1; }
cp -- "$tmp/predictions.out" "$output"
