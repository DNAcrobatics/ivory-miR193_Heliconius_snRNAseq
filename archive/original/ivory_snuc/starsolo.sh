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
