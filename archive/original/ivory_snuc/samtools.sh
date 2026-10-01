#!/bin/bash
#SBATCH -J filterC
#SBATCH -o filterC.out
#SBATCH -e filterC.err
#SBATCH -p nano
#SBATCH -n 32
#SBATCH -t 30:00
#SBATCH --mail-type=END
#SBATCH --mail-user=jjhanly@gwu.edu

##### first segment of code makes the filtered file

#module load samtools
#module load bamtools

#samtools view -@64 -S -b results.Hm_iv_20_Martin_20220914_S29Aligned.out.sam > iv_20.bam
#samtools view -@64 -S -b results.Hm_iv_33_Martin_20220914_S28Aligned.out.sam > iv_33.bam
#samtools view -@64 -S -b results.Hm_iv_38_Martin_20220914_S25Aligned.out.sam > iv_38.bam
#samtools view -@64 -S -b results.Hm_iv_55_Martin_20220914_S26Aligned.out.sam > iv_55.bam
#samtools view -@64 -S -b results.Hm_iv_61_Martin_20220914_S27Aligned.out.sam > iv_61.bam

#list='iv_20 iv_33 iv_38 iv_55 iv_61'
for i in $list
do
#samtools sort -@6 $i.bam -o $i.sorted.bam
#samtools index $i.sorted.bam
#samtools view -@64 -o $i.CortexLocus.bam results.$i*Aligned.sortedByCoord.out.bam Hmel215003o:1400000-1600000
#samtools index $i.CortexLocus.bam
#samtools view -@ 15 -F 0x04 -f 0x2 -q 20 -b $i.CortexLocus.bam > $i.CortexLocus.filt.bam
#samtools view -b -F 16 $i.CortexLocus.bam > $i.CortexLocus.fw.bam
#samtools index $i.CortexLocus.fw.bam
#samtools view -b  -f 16 $i.CortexLocus.bam > $i.CortexLocus.rv.bam
#samtools index $i.CortexLocus.rv.bam
samtools view -h $i.CortexLocus.fw.bam | awk -v OFS="\t" '$0 ~ /^@/{print $0;next;} $6 ~ /N/' | samtools view -b - > $i.CortexLocusSplice.fw.bam
samtools view -h $i.CortexLocus.rv.bam | awk -v OFS="\t" '$0 ~ /^@/{print $0;next;} $6 ~ /N/' | samtools view -b - > $i.CortexLocusSplice.rv.bam

samtools index $i.CortexLocusSplice.rv.bam
samtools index $i.CortexLocusSplice.fw.bam

done
