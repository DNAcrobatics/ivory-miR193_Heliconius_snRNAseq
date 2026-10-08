# Saved workspace export audit — 8 October 2026

The user ran the compact exporter and uploaded [Joe/ivory_workspace_audit](https://drive.google.com/drive/folders/105WFazSqCVKPAH7AW7U1njaQ7-KS1ZQ2). Its inventory lists 225 objects. We inspected selected exported tables and metadata, not the 14 GB binary or expression matrices. Originals and manifests are in `data/workspace_audit/`; command-record text is losslessly gzip-compressed. `metadata/workspace_recovered_files.tsv` records original and stored hashes.

## Recovered saved states

| Object | Nuclei | Saved default assay | State |
|---|---:|---|---|
| Hm_merged | 11,743 | RNA | 17 samples; 6,000 integrated features |
| Hm_sops | 7,035 | SCT | Distinct pre-reintegration lineage branch |
| Hm.sops.sct | 7,035 | integrated | 5,000 integrated features; RNA_by_Ivory present |
| test | 2,649 | RNA | Early scales; identities wt/mut |
| test1 | 973 | RNA | Late scales; identities hi/lo |
| test2 | 1,680 | RNA_by_Ivory | SOP cells; separate normalized-assay branch |

The saved lineage contains 4,397 scale, 1,680 SOP and 958 socket cells. Hm_merged stages are iv_20=32, iv_33=35, iv_38=56, iv_55=46 and iv_61=36. These are recovered values, not confirmation of the manuscript's corrected stages. `metadata/workspace_saved_samples.tsv` preserves them. The broad lineage labels iv_33 as wt in `mut`; its biological heterozygote annotation remains unresolved.

**Do not reuse parent high/low labels for the contrast subsets.** Relative to Hm.sops.sct, test differs at 449 cells, test1 at 102 and test2 at 528. Early scales comprise Hmm_25, Hmr_25, Hmr_30, iv_20 and iv_61; late scales comprise Hmm_40, Hmr_40 and iv_55. Early hi/lo = 1,927/722; late = 832/141. Within early hi, wt/mut = 1,169/758; within late hi = 367/465. Frozen membership files describe these saved states; the DE tables do not record their executed assay, test or cell IDs, so exact invocation is still unproven.

## Saved DE and plotting evidence

| Saved table | Rows | Matching ggplot data |
|---|---:|---|
| WtVMutInHi.Markers | 10,588 | a1 |
| HiVsLo.markers | 7,345 | b1 |
| WtVMutInHiLate.Markers | 11,556 | a4 |
| HiVsLoLate.markers | 10,093 | b4 |

All four match ordered gene IDs and all five statistical columns exactly. `metadata/workspace_plot_matches.tsv` records the checks. HiVsLo.Markers (uppercase M) has the same statistics after gene alignment but different ordering and no plot annotations. WtVMutInHi.Markers2 is absent. Saved UpGenePANGEA has 27 rows and DownGenePANGEA 10, consistent with the late _3 GO/PDF branch already recovered; enrichment generation/background remains missing.

Exact statistics are retained in object_0067, object_0073, object_0221 and object_0222 TSV gzip files. Gene IDs are in `genes`; `.saved_rowname` in annotated tables can be numeric plotting row identifiers.

## Command-record settings and limits

Hm_merged records PCA=50, UMAP/neighbor PCs=1:4, resolution=1, PCA/UMAP seed=42 and clustering seed=0; its integrated assay has 6,000 features. Hm.sops.sct has 5,000 features and records PCA=50, UMAP PCs=1:3, latest neighbors PCs=1,2,3,4,6, resolution=1 and the same recorded seeds. Its clustering timestamp precedes the latest neighbor and UMAP timestamps: the latest neighbor command cannot be assumed to have generated its saved clusters. These records do not establish the manuscript's 3,000-feature/12-PC/resolution-0.5 branch. Existing historical reconstruction scripts are not silently changed to a competing branch.

The export session reports R 4.5.1, Seurat 5.5.0, SeuratObject 5.4.0 and sctransform 0.4.3 on macOS. This is the export environment, not necessarily the original analysis environment.

## Restore the frozen labels

`scripts/snrna/07_restore_saved_metadata.R` selects saved full cell IDs and restores biological fields without recalculating thresholds. It requires explicit assay and identity choices, rejects absent cells/library mismatches, and cannot reconstruct a missing RNA_by_Ivory assay.

```bash
Rscript scripts/snrna/07_restore_saved_metadata.R \
  --input inputs/Hm.sops.sct.rds \
  --metadata data/workspace_audit/object_0196_metadata.tsv.gz \
  --assay RNA --identity mut --output results/early_saved_state.rds
```

Files `early_hi_vs_lo.membership.tsv`, `late_hi_vs_lo.membership.tsv`, `early_wt_vs_mut_in_hi.membership.tsv` and `late_wt_vs_mut_in_hi.membership.tsv` can feed 05_find_markers.R once actual DE assay/test/filter settings are established. Metadata alone does not establish these settings.

## Remaining completion inputs

1. Confirm final panel/branch selection and DE execution settings: assay, test/filter parameters, group direction and nucleus versus replicate-level analysis. A successful call/log is useful; the four recovered tables do not need resending.
2. Confirm corrected stages and iv_33 treatment. Saved per-cell metadata and identities are recovered; no repeat export is requested.
3. Supply expression-object/count-matrix provenance and upstream SCT/QC/import setup for scientific reruns. Metadata cannot recover expression values or omitted history. RNA_by_Ivory construction and the final classification rule still need provenance.
4. Reference/index/BAM commands and locus GFF, RNAhybrid inputs/full targets, collaborator attribution, enrichment background/settings, license and release DOI remain as detailed in FILE_REQUESTS.md. Supplementary tables and figure-assembly locations are now recovered; see PAPER_DOCUMENTS.md.

## Verification

Downloaded files covered by MD5SUMS.tsv matched every listed checksum. Metadata cell IDs were unique; early/late cohorts are contained in the saved lineage; four DE tables matched saved plotting statistics exactly. SHA256 hashes were recorded. Rscript is unavailable here: the new recovery script was reviewed but not executed. These data checks do not validate the scientific model or final publication mapping.
