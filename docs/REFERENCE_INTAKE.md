# Reference intake and narrowed requests — 8 October 2026

**Paper-document scope update — 8 October 2026:** The manuscript/reporting documents, all eight supplementary tables and figure-assembly locations are recovered from Drive; no repeat upload is requested. Nanopore was collaborator work and is outside the user code-intake request. See [PAPER_DOCUMENTS.md](PAPER_DOCUMENTS.md). S8 supplies the target/overlap annotation table and 15 conserved candidates; remaining miR requests concern prediction inputs and preparation/parsing code.

The user supplied the ribosomal and mitochondrial gene lists, the reduced-Ivory counting GTF and Luca-reannotated locus GFF. Originals are preserved byte-for-byte, with the large GTF losslessly compressed in `data/references/`. Original and stored SHA256 checksums are in `metadata/reference_intake_20261008.tsv`.

## Reference checks and use

The ribosomal list has 160 rows representing 150 unique IDs; the mitochondrial list has 38 unique IDs. All 150 ribosomal IDs and all 38 mitochondrial IDs occur as exon gene IDs in the supplied GTF. Preserve the duplicate ribosomal rows as source evidence; PercentageFeatureSet uses matching feature IDs, not row multiplicity. These are annotation matches, not validation against every expression object.

The GTF has 506,772 nine-column records (262,541 exon and 244,231 CDS records), with 18,130 distinct exon gene IDs. Every checked interval had positive, ordered coordinates and every exon had a gene_id attribute. This structural check does not establish reference FASTA compatibility or reproduce STAR alignment.

**The two Ivory annotations serve different purposes.** The counting GTF defines one Ivory exon at Hmel215003o:1430470–1432810, whereas the locus GFF gene spans 1430436–1568488 and includes transcript/exon structure. The locus GFF has 1,180 nine-column records with valid intervals. It contains repeated IDs, including gene:Ivory used for both gene and mRNA; it is preserved unchanged and is not asserted to pass strict GFF3 validation. Do not replace the counting annotation with this locus visualization annotation.

Decompress the counting GTF before passing it to STAR:

```bash
gzip -dc data/references/Hmel.v3.1.annotation.CAT.MOD.plusMito.reducedIvory.gtf.gz > inputs/Hmel.v3.1.annotation.CAT.MOD.plusMito.reducedIvory.gtf
```

Create `inputs/` first and choose a new destination if a file already exists. QC scripts can use `--ribo data/references/Hmel_ribosomal_genes.txt --mito data/references/Hmel_mito_genenames.txt`. The locus GFF is now available for reconstructing the gene-model panel; the current read-start plotter does not yet draw annotation tracks. No R, STAR or BAM processing was rerun here.

## Existing reanalysis source: supplied, not requested again

The user identifies `archive/drive/iv_reanalysis/ivory_reanalysis(1).r` as the analysis and metadata source. It supplies lineage selection, SCT/reintegration trials, biological labels and DE candidate calls. The missing-file request for this code is closed; remaining questions concern which executed branch produced the final outputs.

The script begins by loading Hm_combined.sct.h5seurat (line 9), rather than importing matrices/building that object. Its lineage trial selects 7,000 integration features (line 147), whereas the saved Hm.sops.sct assay has 5,000; both are retained as distinct evidence. Its mutation labels select iv_61, iv_38 and iv_55 as mut and label the remainder wt (lines 235–239), explaining the saved iv_33 label. Biological heterozygosity and final plot treatment must still be distinguished.

The early WtVMutInHi.Markers call at line 1010 operates on all of test, without a hi subset despite the variable name. Its high/low call follows at line 1036. Other supplied versions/history use different subsets/calls. Therefore do not infer high-only membership from a table name or silently rewrite the original executed comparison. Frozen memberships from the saved export remain candidate inputs pending result provenance. Original source and cleaned workflows are already supplied; a new duplicate script upload is unnecessary.

## Narrow miR upload request

The listing establishes the main directory as `/gpfs/automountdir/gpfs/scratch/martinlab/joe/miR`. File contents have not yet been received. The smallest useful bundle is:

- `rnahybrid.sh`
- `hme-miR-193.fa`
- `HmelFrancGenes.longest.3utr.fa`
- `longest_transcript_ids.txt`
- `Hmel_miR193_RNAhybrid_longest3utr.out`
- `Hmel_miR193_RNAhybrid_gene_hit_summary.tsv`
- `Hmel_miR193_RNAhybrid_predicted_genes.txt`
- `rnahybrid_mir193.*.out` and `.err`
- Any scripts/history that selected the longest transcripts, extracted 3′ UTRs and converted predictions to the gene summary/with_miRhit.

The final UTR FASTA is only 8.9 MB, so the whole 640 MB directory and intermediate 37–80 MB FASTAs are unnecessary for this first intake. First_pass_with_miTAR is needed only if that branch contributed to the paper or contains the preparation/parsing scripts. Cross-species comparison inputs and Nanopore files remain separate requests.
