# Remaining files for the ivory mir-193 paper code repository

Updated 2026-10-01 after inspecting the manuscript, supplementary figures, reporting summaries, editorial checklist, two pasted analysis versions and the inline source. The broad first-pass requests have now been narrowed to the files referenced by the recovered code. No need to resend the manuscript or `CellCycleGenes.txt` already supplied.

## Received and incorporated

| Material | What we now have | Remaining limitation |
|---|---|---|
| `FinalEdits_Ivory_d78k_scRNAseq_.docx` | Current title, main figure legends, Methods, data/code statements and target-prediction Methods | Contains internal discrepancies listed below; not altered |
| `Supplementary_dta.docx` | Supplementary figure legends and embedded figures S1–S8 | This file does not contain the machine-readable DE/enrichment/target tables S4–S8 |
| Author checklist and both reporting PDFs | Editor requirements, manuscript-reported software versions and BioProjects | Runtime version logs and a specific repository/release DOI still needed |
| `ivory_reanalysis.r` terminal paste | Original body extracted with shell prompts removed; exact transcript retained | Exploratory source with unresolved expressions/session dependencies |
| Longer standalone pasted analysis | Figure 5–7, S7/S8 candidates, GO plotting, high/low normalization attempt and miRNA highlight lists | Original filename unconfirmed; multiple candidate versions; no inputs/results supplied |
| `starsolo.sh` | STARsolo command, read order, barcode/UMI lengths and counting settings | `jobs`, edited GTF, reference build, logs and final BAM-producing command missing |
| `samtools.sh` | Locus/strand commands in comments plus active spliced-read extraction | `list` assignment and prerequisite commands commented out; not a complete locus-track workflow |
| `wt_vs_ivory_new.May2024.r` | Historical import/QC/integration for ten samples | Creates `Iv_combined.full.sct`, not the `Hm_combined.sct.h5seurat` loaded by later code |
| `CellCycleGenes.txt` | Full 66-row mapping table, including original `NO_HIT` and duplicated Hmel IDs | Need actual feature IDs used/returned warnings; source table preserved unchanged |

Original snapshots and their checksums are in `archive/` and `metadata/received_source_manifest.tsv`. They are not yet executable publication workflows. Input/output filename references have been indexed separately.

## Send next — highest priority

### 1. Code that created the combined object

Send the script(s) that wrote **`Hm_combined.sct.h5seurat`**, including full wild-type time-series import, updated Piano Keys metadata, QC and integration. This is the main upstream dependency for both pasted analysis versions. It may be in `ivory_proj_transfer_aug2025`, `iv_v_wt_timeseries_March2025` or `iv_relignment_reanalysis`; search for the exact object filename. Also send any script that generated a final **`Hm.sops.sct.h5seurat`**, `Hm_sops.h5seurat` or an equivalent RDS used for the paper. Initially, an object filename/size manifest plus its per-cell metadata is enough; large objects can follow if needed.

Send these already-listed source files from `Iv_v_Wt_comp` as well:

- `rscript.loop.r`
- `rscript.r`
- `wt_vs_ivory_scales.r`

The May 2024 script is already transcribed. Supplying actual script files in a bundle later will let us verify transcription and recover timestamps, but is not the next blocker.

### 2. Small input files explicitly referenced by the figure code

From `ivory_proj_transfer_aug2025` (or wherever the script was run):

- `Hi_wt_v_mut_for_dotplot.txt`
- `WtVMutInHiLate.Markers_genenames.txt`
- **`HiVsLoLate.markers_naems.txt`** — retain this spelling; it is what the script reads
- `HomologGeneNames.txt`
- `fbgn_fbtr_fbpp_expanded_fb_2024_01.tsv`
- `UpGenes.txt`, `DownGenes.txt`, `UpGenesLate.txt`, `DownGenesLate.txt`
- `UpGenes.txt_FBpp.list.txt`, `DownGenes.txt_FBpp.list.txt`
- `UpGenesLate.txt_FBpp.list.txt`, `DownGenesLate.txt_FBpp.list.txt`
- `UpGenes_DirectGoEnrichment_PANGEA.csv`
- `DownGenes_DirectGoEnrichment_PANGEA.csv`
- `UpGenesLate_DirectGoEnrichment_PANGEA_3.csv`
- `DownGenesLate_DirectGoEnrichment_PANGEA_3.csv`
- `figure_geneCodes.txt` if its plotting branch was retained
- `gene_names_for_volcplot.txt` only if the commented older branch was used

From the older import/reference folders:

- `Iv_v_Wt_comp/Hmel_ribosomal_genes.txt`
- `Hmel2.5_PLUS_MITO_GENOME/Hmel_mito_genenames.txt`
- `ivory_snuc/jobs`
- `Hmel2.5_PLUS_MITO_GENOME/Hmel.v3.1.annotation.CAT.MOD.plusMito.reducedIvory.gtf`
- The script that generated that edited GTF and the STAR index; the reference FASTA filename, source/version and checksum

The gene-set files must include the selection/export rule used to create them. The pasted code does not establish a reliable final significance threshold for all four enrichment inputs.

### 3. Actual final statistical and plotting outputs

Send the machine-readable Supplementary/Source Data tables **S4–S8** and the complete saved DE tables for the four plotted comparisons. The pasted R variable names are:

- Early high vs low: `HiVsLo.markers`
- Late high vs low: `HiVsLoLate.markers`
- Early WT vs mutant within high ivory: `WtVMutInHi.Markers` (also a newly computed `WtVMutInHi.Markers2` appears)
- Late WT vs mutant within high ivory: `WtVMutInHiLate.Markers`

If these only exist in an R session, export them with rownames/gene IDs preserved. Include all columns, notably `p_val`, `p_val_adj`, `avg_log2FC`, `pct.1` and `pct.2`. Include the actual per-cell metadata and identities used for each comparison, not only the gene table. A frozen plot-input table can distinguish the final output from nearby exploratory comparisons without rerunning everything.

Final panel PDFs from the recovered code would help verify the mapping:

- `Hm.sops.sct.celltypes.pdf`, `Hm.sops.sct.race.pdf`, `Hm.sops.sct.time.pdf`
- `RaceClusters.pdf`, `PhaseClustersEarly.pdf`, `PhaseClustersLate.pdf`
- `EarlyScaleCluster_IvoryFeature_and_origIdents.pdf`, `LateScaleCluster_IvoryFeature_and_origIdents.pdf`, `All_IvoryFeature.pdf`
- `VolcPlot_WtVMutInHi.pdf`, `VolcPlot_WtVMutInHi_late.pdf`
- `DotPlot_WtVMutInHi.pdf`, `DotPlot_WtVMutInHi.Late.pdf`
- `DotPlot_HiVsLo.pdf`, `DotPlot_HiVsLo.Late.pdf`, `VolcPlotBoth.pdf`, `VolcPlotBothLate.pdf`
- `GoTermsOut.pdf`, `GoTermsOutLate.pdf`
- `DotPlot_allHm_merged.pdf`, `FeaturePlots_allHm_merged.pdf`, `DimPlots_allHm_merged.pdf`
- `Fig_Sx_gene_counts_in_scale_cells_histogram.pdf`

Use the final versions even if renamed. `VolcPlot_HiVsLo.pdf` is reused for early and late output in the pasted script, so that basename alone cannot identify both panels.

## Other required branches

### Locus readthrough and Figure 6A

Send the **cell-barcode selection/filtering script**, **forward-strand coverage/depth script**, **ribbon-track R script**, relevant locus BED/GTF models and the small coverage tables used for the final panel. The provided `samtools.sh` selects CIGARs containing `N`; this is spliced-read extraction, not proof of the full exonic/intronic coverage pipeline. The STARsolo command also does not contain the active BAM/CB output settings needed to reconstruct these tracks. Include the actual BAM-producing command, barcode-tag handling, track units/normalization and final per-sample cell membership.

### RNAhybrid and comparative target analysis — confirmed in this manuscript

Send the launch and parsing/collapse scripts, mature miR-193 FASTA, **`HmelFrancGenes.longest.3utr.fa`** if that is still the actual target FASTA (earlier filename lead), its generation code and reference manifest, full prediction output and the final target table. The new manuscript specifies RNAhybrid **2.1.2**, **`-s 3utr_fly`** and **`-e -20`**, but the original command and sequence inputs are not supplied yet.

Send the code/data that create **`with_miRhit`**, the provenance of **`top15_mir_multi`**, the gene-ID/orthology mappings, and the source/version of the **Bicyclus** and **Junonia** DE tables. The hardcoded 15-gene list is present; its derivation is missing. Include the overlap selection thresholds and direction rules.

### Enrichment generation and correction

The recovered code **reads PANGEA-named DirectGO result CSVs** and plots them; it does not generate those results. Send the original PANGEA/DirectGO job exports/settings or generating code, input lists, gene universe/background, annotation/version and full unfiltered result CSVs. Supply modPhEA exports/settings if used in final tables or Results. The figure legend mentions permutation-based significance, but the pasted plotting code alone does not verify the underlying test. Several plotting blocks mix **Benjamini–Yekutieli** filtering with **Benjamini–Hochberg** coloring; final CSVs and figure inputs are needed to establish which adjustment belongs to each panel.

### Nanopore crispant genotyping and Figure S6

Send Dorado basecalling/demultiplexing commands with the full model name, barcode/sample map, minimap2/samtools commands, the custom R coverage plotter and its depth tables. `ivory_snuc/samtools.sh` concerns snRNA-seq locus reads; it does not supply the Nanopore genotyping pipeline. Primer/guide and genotype summary tables S1–S3 should accompany the machine-readable supplement.

### Metadata, environments and release details

- Final sample/library/run sheet: individual vs library IDs, genotype, stage, population, exclusions and SRA run accessions. Preserve old sample IDs such as `iv_61` even when the stage has been corrected.
- Per-cell metadata with barcode, `orig.ident`, stage, race/population, genotype, cell type, high/low ivory, phase, final assay identity and contrast membership. Include final cell counts per sample/group.
- Actual final R `sessionInfo()`, Seurat/SeuratObject/SeuratDisk/sctransform versions, environment exports and STAR/samtools logs. The reporting forms say STAR **2.7.11b**, samtools **1.2.1**, Seurat **5.4.0**; these are reported values, not runtime checks.
- Helper startup/functions/packages defining `%ni%`, `scale_color_batlow`, `scale_colour_discreterainbow`, `colours` and any palettes used in the final plotting session.
- One successful run log per branch, seeds and expected output summaries.
- Desired repository name under **Hanliconius**, code contributor credits and code license choice. The manuscript currently points to the GitHub profile rather than a specific repository.
- A tagged code release archived in a DOI-minting repository. Your supplied editorial checklist explicitly requests this, so GitHub plus an archived release should be the final delivery plan.

## Issues to settle with the final source and outputs

1. **DE replication unit:** the main figure candidates use `FindMarkers` on nuclei. The `AggregateExpression` blocks supplied are for expression ranking/histograms, not a replicate-level DE model. The legends specify Wilcoxon/Bonferroni while Results/Methods say pseudobulk. Preserve the actual analysis and make those descriptions consistent after its final version is identified.
2. **High/low classification and assay:** one branch uses RNA thresholds of `Ivory > 1` for Piano Keys and `> 4` for other populations. Another creates `RNA_by_Ivory` using per-sample mean ivory expression, then ends with an incomplete threshold expression. Need the final assignment and which assay supplied DE; do not choose between these branches from prose alone.
3. **Metadata and integration versions:** the May 2024 code uses stages 32/35/58/46/36 and 5,000 features/8 PCs/resolution 1. The new manuscript uses revised stages and 3,000 features/12 PCs/resolution 0.5; the lineage paste also contains 7,000 features and 3 PCs. Request the exact final object-building code and metadata.
4. **Figure 7 cell fate label:** its introductory legend says high-ivory cells acquire Type I in wild type, whereas Results describe presumptive Type II cells converting to Type I in mutants. This affects contrast documentation; verify the intended wording against the final memberships.
5. **Genotype labeling:** the May script includes `iv_33` as a heterozygote, while an exploratory later block labels every sample not in three homozygous mutants as `wt`. Main early/late subsets exclude `iv_33`, but broader plots need explicit genotype membership.
6. **Source completeness:** missing parentheses, bare notes, shell commands inside R and session-only objects prevent sourcing the long paste directly. Originals are retained; cleanup will follow input recovery, with changes logged and results compared.

## Efficient collection command

Download `collect_project_code.py` from the earlier bundle and run from its directory. Edit folder names if different. This command includes text/CSV/TSV lookup and plotting tables as well as code, up to 20 MiB per file; larger files are inventoried for later retrieval. It does not run analyses or change source files and refuses to overwrite an existing bundle.

```bash
python3 collect_project_code.py \
  --root /gpfs/automountdir/gpfs/scratch/martinlab/joe/ivory_proj_transfer_aug2025 \
  --root /gpfs/automountdir/gpfs/scratch/martinlab/joe/ivory_snuc \
  --root /gpfs/automountdir/gpfs/scratch/martinlab/joe/Iv_v_Wt_comp \
  --root /gpfs/automountdir/gpfs/scratch/martinlab/joe/iv_relignment_reanalysis \
  --root /gpfs/automountdir/gpfs/scratch/martinlab/joe/iv_v_wt_timeseries_March2025 \
  --root /gpfs/automountdir/gpfs/scratch/martinlab/joe/Hmel2.5_PLUS_MITO_GENOME \
  --include-glob '*.txt' \
  --include-glob '*.csv' \
  --include-glob '*.tsv' \
  --include-glob '*.gtf' \
  --out "$PWD/ivory_paper_code_collection_02.tar.gz"
```

Add a root for the RNAhybrid and Nanopore work if those are elsewhere. `he_ts` remains conditional on actual use in this paper. Large raw reads, BAMs and Seurat objects are not the first upload request; send filenames/sizes/accessions and per-cell metadata first.

## Original request IDs and current status

| ID | Status after this intake | Still required |
|---|---|---|
| R01 Manuscript/figures | Manuscript, main legends and supplement S1–S8 received | Source-data tables, final exported panel files/assemblies, repository/release details |
| R02 Code | Five analysis/script snapshots received | Upstream object builder, remaining candidates, readthrough, RNAhybrid and Nanopore scripts |
| R03 Samples | Historical May 2024 values recovered | Final corrected sample/run and per-cell metadata |
| R04 References | Exact reduced-ivory GTF path and mapping dependencies identified | Edited GTF, build code, FASTA manifest, gene/orthology tables |
| R05 Alignment | STARsolo snapshot received | `jobs`, final reference/index build, version/logs and BAM/barcode output command |
| R06 Seurat/cell cycle | Lineage/plotting code and cell-cycle mapping received | Final combined-object builder, final assignments/objects, helper definitions |
| R07 DE | Four main contrast candidate calls identified | Final assay/threshold/replicate description, saved full tables and group metadata |
| R08 Enrichment | Plotting code and four CSV basenames identified | Actual CSVs and their generation/settings/background/correction |
| R09 Locus tracks | Partial historical samtools commands received | Barcode selection, full depth/strand pipeline and final coverage/ribbon inputs |
| R10 Targets | Manuscript Methods and 15-gene highlight list received | RNAhybrid inputs/code, `with_miRhit`, cross-species source tables/derivation |
| R11 Nanopore | Manuscript description received | Actual processing/plotting scripts and depth/sample tables |
| R12 Figure exports | Many exporter names recovered | Final inputs/PDFs, manual assembly steps and mapping confirmation |
| R13 Environments | Reported versions recovered | Actual session/package/tool records and environments |
| R14 Validation | Static source inspection only | Successful logs and expected outputs; no scientific rerun performed |
| R15 Data/release | BioProjects PRJNA1404329 and PRJNA1118807 provided | Run-level mappings, processed-input locations, code license/release archive DOI |
