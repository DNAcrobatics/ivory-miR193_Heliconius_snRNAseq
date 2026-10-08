#!/bin/bash
#SBATCH --job-name=rnahybrid_mir193
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=rnahybrid_mir193.%j.out
#SBATCH --error=rnahybrid_mir193.%j.err

micromamba activate rnahybrid

RNAhybrid \
  -t HmelFrancGenes.longest.3utr.fa \
  -q hme-miR-193.fa \
  -s 3utr_fly \
  -m 10000 \
  -b 10 \
  -e -25 \
  > Hmel_miR193_RNAhybrid_longest3utr.out
