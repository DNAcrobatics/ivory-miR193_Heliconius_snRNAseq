# Saved-workspace audit — 6 October 2026 (New York)

**8 October 2026 update:** The compact export has been received and inspected. [SAVED_WORKSPACE_AUDIT.md](SAVED_WORKSPACE_AUDIT.md) documents recovered metadata, assay summaries, four plotted DE tables and frozen candidate memberships. Earlier requests below for another compact export, these tables or saved per-cell metadata are superseded. Final execution provenance and scientific reruns remain pending.

## Upload and inspection status

The refreshed Drive listing exposes [.RData](https://drive.google.com/file/d/1FKwR1xzCkHQ4AIl9E9lU7SXYYu_Q3oPh/view?usp=drivesdk) at **15,360,101,380 bytes** (14.31 GiB). This agrees with the rounded 14G local listing. The connector supplied neither an upload-progress field nor the requested MD5 checksum; the full binary has not been downloaded, checksum-validated, decompressed or loaded. Its presence and size are confirmed, not its loadability or saved object contents.

`.RDataTmp` is listed at 3,989,831,680 bytes, substantially smaller than its local 7.6G listing; treat it as a potentially incomplete failed copy. `.RDataTmp1` is listed at 8,132,533,918 bytes. Neither temporary file was used for reconstruction.

The refreshed `.Rhistory` is 29,407 bytes and its SHA-256 is `37907d297559266c7c474fbd712f5d7f4b41464ef92025377f0faae4bfb98b74`, identical to the previously archived copy. It is a command record, not an execution log: it includes failed expressions and does not establish which commands succeeded.

This inspection runtime has **no R/Rscript executable**, an **8 GiB cgroup memory limit**, and no swap. The compressed workspace is already larger than that memory allowance. A full Seurat workspace load would not be a viable audit here; no object counts, final assay identities, memberships or DE tables are claimed as recovered from it. Disk space alone would not solve this limit.

## Findings from history and saved GO output

Line numbers refer to `archive/drive/iv_reanalysis/Rhistory.txt`; the evidence index is `metadata/history_audit_evidence.tsv`.

- Lines 1–146 contain repeated SOP UMAP/neighbor/cluster experiments. Lines 149–165 include 5,000 SCT variable features, time regression, five PCs for UMAP/neighbors, resolution 2 and package reinstall attempts. These are additional historical trials, not proof of final manuscript settings or successful package installation.
- Lines 310–316 record RNA high/low thresholds of `Ivory > 1` in race `iv` and `> 4` elsewhere, consistent with the existing cleaned classification branch. This does not settle the later normalized `RNA_by_Ivory` experiment or prove the final objects use these labels.
- Line 319 sets `Idents(test1)` to the high/low `Iv` labels, then line 321 asks `FindMarkers(..., ident.1='wt')`. There is no intervening genotype identity assignment. Lines 337 and later repeats have the same unresolved state. This sequence cannot establish a successful WT-versus-mutant comparison; the saved tables and identities are needed. The cleaned DE script already uses two explicit groups and frozen membership instead.
- Lines 298, 375 and 399 record RNA high-versus-low FindMarkers calls with `recorrect_umi=FALSE` and omitted `test.use`. The history gives no replicated pseudobulk DE implementation. The statistical description remains unresolved until the executed workflow/version and tables are recovered.
- GO exploration alternates Bonferroni, BH and BY filters. The `_2` inputs occur in earlier trials. At lines 433–508 both final panels are rebuilt from `_3`: BY < 0.05 filtering for both, BY coloring for upregulated terms and BH coloring for downregulated terms. This agrees with the larger script's `_3`/`_3` candidate.

The accessible [GoTermsOutLate.pdf](https://drive.google.com/file/d/1RsOU4yHURHTGi8BHqZ9IDRwuzjoWNE6D/view?usp=drivesdk) was downloaded, text-extracted and visually reviewed. Its top panel has **23 distinct term labels**, consistent with the **27 BY-significant rows** in `UpGenesLate_DirectGoEnrichment_PANGEA_3.csv` (some term labels recur). Its bottom panel has **10 term labels**, consistent with `DownGenesLate_DirectGoEnrichment_PANGEA_3.csv` under BY < 0.05. The `_2` files supply only 3 up and 6 down labels and do not explain the complete saved panel. The history and saved PDF therefore both support the `_3`/`_3` candidate for this particular output. It is a term-label comparison, not a numerical or final assembled-paper validation. Preserve all 27 up rows rather than deduplicating by display label.

`metadata/go_export_audit.tsv` records row counts, correction-specific counts and checksums across all eight recovered CSVs. Filter/color correction choices remain explicit in the cleaned plotter.

## Smallest useful export

Use `tools/export_workspace_audit.R` in an **existing R session that already has the workspace loaded**. This avoids creating a second full copy in RAM. Download that helper from the repository, then run in R:

```r
source("export_workspace_audit.R")
export_workspace_audit(outdir = "ivory_workspace_audit")
```

If the helper was obtained with the repository structure intact, source `tools/export_workspace_audit.R` instead. Upload the resulting `ivory_workspace_audit` folder into the shared Drive folder. The export contains an object inventory, saved table objects, per-cell metadata and active identities, assay/layer summaries, recorded Seurat command parameters, saved ggplot input and small gene/target vectors. No count matrices or full Seurat objects are exported; tables are gzip-compressed and large tables are flagged rather than copied without limit. `exports.tsv` records skipped files and export errors. It does not recompute science or infer contrast groups.

For a fresh process on a sufficiently provisioned machine with compatible R/Seurat packages:

```bash
Rscript --vanilla tools/inspect_workspace.R /path/to/.RData ivory_workspace_audit
```

This still loads the full workspace once and needs enough RAM. The compressed file size is not a memory estimate. The source on the Mac was the `.RData` in the terminal's `iv_reanalysis` working directory; its absolute source path was not supplied. The Drive-mounted destination is not a substitute for that original local file when disk space is tight.

The export helper was manually reviewed here but could not be executed without R. The inspection sessionInfo describes the export environment, not necessarily the original analysis environment. Saved active identities and recorded commands describe saved state; they do not prove earlier identities used by each DE call.

## Completion status

The workspace link is recovered; do not resend it. Next priority is this compact export, especially `test1`, `Hm.sops.sct`, `Hm_merged`, `WtVMutInHiLate.Markers`, `HiVsLoLate.markers` and any early equivalent result tables. Other outstanding branches remain reference construction/count inputs, RNAhybrid targets, Nanopore processing, final figure mapping, license and release archive. Scientific reruns and final result matching remain pending.
