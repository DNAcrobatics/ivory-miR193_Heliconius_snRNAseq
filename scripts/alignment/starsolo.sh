#!/usr/bin/env bash
# Clean starsolo.sh: same counting settings; explicit paths and sample selection.
# Usage: bash starsolo.sh jobs INDEX fastq_dir genome_dir annotation.gtf out_dir THREADS
# INDEX is zero-based; on Slurm pass "$SLURM_ARRAY_TASK_ID" explicitly.
# Activate the STAR environment before calling; scheduler options belong in sbatch.
set -euo pipefail
fail() { printf '%s\n' "$*" >&2; return 1; }
[[ $# == 7 ]] || fail 'Usage: starsolo.sh jobs INDEX fastq_dir genome_dir annotation.gtf out_dir THREADS'
jobs=$1; task_index=$2; fastq_dir=$3; genome_dir=$4; annotation=$5; out_dir=$6; threads=$7
[[ $task_index =~ ^[0-9]+$ && $threads =~ ^[1-9][0-9]*$ ]] || fail 'INDEX must be >=0 and THREADS >0'
task_index=$((10#$task_index))
[[ -f $jobs && -f $annotation && -d $genome_dir ]] || fail 'Missing jobs/reference/annotation'
command -v STAR >/dev/null || fail 'STAR is not available in the active environment'
# Ignore blank lines, but reject names with whitespace or path components.
mapfile -t names < <(awk 'NF {print}' "$jobs")
[[ $task_index -lt ${#names[@]} ]] || fail 'INDEX exceeds jobs list'
sample=${names[$task_index]}
[[ $sample =~ ^[A-Za-z0-9_.-]+$ ]] || fail 'Invalid sample name'
r2="$fastq_dir/${sample}_R2_001.fastq.gz"; r1="$fastq_dir/${sample}_R1_001.fastq.gz"
[[ -f $r2 && -f $r1 ]] || fail 'Missing FASTQs'
mkdir -p "$out_dir"
prefix="$out_dir/redo.$sample."
[[ ! -e ${prefix}Log.out && ! -d ${prefix}Solo.out ]] || fail 'STAR output prefix already exists'
printf 'Sample: %s\n' "$sample"
STAR --runThreadN "$threads" \
  --genomeDir "$genome_dir" --outFileNamePrefix "$prefix" \
  --readFilesCommand zcat --readFilesIn "$r2" "$r1" \
  --alignIntronMax 75000 --sjdbGTFfile "$annotation" \
  --soloType CB_UMI_Simple --soloCellFilter EmptyDrops_CR \
  --soloFeatures GeneFull --soloMultiMappers EM --soloCBwhitelist None \
  --soloUMIlen 12 --soloCBlen 16 --quantMode GeneCounts \
  --soloBarcodeReadLength 0 --sjdbGTFtagExonParentGene gene_id
# Historical BAM/barcode options were commented; they are not silently enabled here.
