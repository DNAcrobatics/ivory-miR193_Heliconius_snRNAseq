# Strand-specific ivory locus coverage

**Paper:** Figure 6A.

**Status:** Partial historical BAM processing.

Recover per-barcode exonic/intronic coverage and the locus ribbon plot.

## Inputs

- Aligned BAMs with retained nucleus barcode tags and indices.
- Final plotted barcode selections, strand convention, annotation and plotting coordinates.
- Existing strand-separated CortexLocus BAMs for the historical script below.

## Outputs

- Historical strand-specific spliced BAMs; final per-cell coverage tables and Figure 6A vector plot still needed.

## Software

samtools, awk and original coverage/plotting software; actual versions pending.

## Execution and interpretation

The historical list assignment is commented, so its loop does nothing unless list is already set in the shell. Earlier extraction, filtering and strand splitting are commented as well. The active commands retain reads whose CIGAR contains N from pre-existing strand-specific BAMs. This is only part of the coverage workflow, and it does not distinguish all exonic/intronic coverage or make the final ribbon plot.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Code

**Source:** [samtools.sh](../archive/original/ivory_snuc/samtools.sh).

**Archive treatment:** supplied snapshot reproduced unchanged. Paths and scheduler settings remain historical; this is not a portable final workflow.

```bash
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
```


[Return to analysis index](../README.md#analysis-index).
