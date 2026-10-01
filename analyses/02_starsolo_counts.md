# Single-nucleus alignment and counting

**Paper:** Upstream of Figures 5–7 and S7–S8.

**Status:** Historical script.

Align snRNA-seq reads and generate STARsolo gene-full count matrices.

## Inputs

- Paired gzipped FASTQs; historical script passes R2 first and R1 second.
- The ordered `jobs` sample list; it has not been supplied.
- The reference from workflow 01.

## Outputs

- STARsolo GeneFull matrices and alignment logs under each `redo.<sample>` prefix.

## Software

STAR, samtools, micromamba and Slurm. Reported versions need confirmation against run logs; see [environments](../envs/README.md).

## Execution and interpretation

The source is a five-task Slurm array. Confirm which run supplied the final filtered matrices and how `redo.*` outputs became `iv_*.filtered`/`filtered.Hm*`. Preserve CB length 16, UMI length 12, no whitelist, EM multimappers, and EmptyDrops_CR unless the final run records say otherwise. BAM and barcode attribute options are commented out, so downstream per-barcode coverage cannot yet be reproduced from this script alone.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Code

**Source:** [starsolo.sh](../archive/original/ivory_snuc/starsolo.sh).

**Archive treatment:** supplied snapshot reproduced unchanged. Paths and scheduler settings remain historical; this is not a portable final workflow.

```bash
#!/bin/bash
#SBATCH -J starsolo
#SBATCH -o starsolo.out
#SBATCH -e starsolo.err
#SBATCH -p nano
#SBATCH -n 64
#SBATCH -t 30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jjhanly@gwu.edu
#SBATCH --array=0-4%5

module load samtools
micromamba activate STAR


names=($(cat jobs))

echo ${names[${SLURM_ARRAY_TASK_ID}]}


mkdir test.results
STAR --runThreadN 64 \
 --genomeDir ../Hmel2.5_PLUS_MITO_GENOME \
 --outFileNamePrefix redo.${names[${SLURM_ARRAY_TASK_ID}]} \
 --readFilesCommand zcat \
 --readFilesIn ../../221004_A00201R_0751_AHCMFLDSX5/fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R2_001.fastq.gz ../../221004_A00201R_0751_AHCMFLDSX5/fastq_files/${names[${SLURM_ARRAY_TASK_ID}]}_R1_001.fastq.gz\
 --alignIntronMax 75000 \
 --sjdbGTFfile ../Hmel2.5_PLUS_MITO_GENOME/Hmel.v3.1.annotation.CAT.MOD.plusMito.reducedIvory.gtf \
 --soloType CB_UMI_Simple \
 --soloCellFilter  EmptyDrops_CR \
 --soloFeatures GeneFull \
 --soloMultiMappers EM \
 --soloCBwhitelist None \
 --soloUMIlen 12 \
 --soloCBlen 16 \
 --quantMode GeneCounts \
 --soloBarcodeReadLength 0 \
 --sjdbGTFtagExonParentGene gene_id \



## removed to make sam
# --outBAMcompression -1 \
# --outSAMtype BAM SortedByCoordinate
#  --outSAMattributes CB


#removed to see if it makes the counts work
# --sjdbGTFtagExonParentTranscript Parent
# --sjdbGTFtagExonParentGene gene \
# --sjdbGTFfeatureExon exon \

# --outSAMtype BAM SortedByCoordinate \
# --outSAMattributes GX \

# --sjdbGTFtagExonParentGene Parent \

# --outSAMattributes CB CR CY UB UR UY \
# --outSAMtype BAM SortedByCoordinate \
# --limitBAMsortRAM 5572187465
```


[Return to analysis index](../README.md#analysis-index).
