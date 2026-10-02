# Count import, QC and integration

**Paper:** Full dataset underlying Figures 5–7 and S7.

**Status:** cleaned scripts supplied; final inputs/settings and result matching pending.

Document matrix import and recover the full integrated object used for the paper.

## Inputs

- Ten historical filtered 10X matrix directories named in the code.
- Ribosomal and mitochondrial gene lists.
- Final sample metadata, all additional time-series matrices, and the builder for `Hm_combined.sct.h5seurat`.

## Outputs

- Historical `Iv_combined.full.sct` object; final `Hm_combined.sct.h5seurat` integration/save tail is now recovered; upstream setup and final-run confirmation remain missing.

## Software

R, Seurat, SeuratDisk, sctransform, ggplot2, dplyr and patchwork; exact package environment pending.

## Execution and interpretation

This May 2024 branch imports ten samples and sets historical times/genotypes. It uses 5,000 variable features, eight PCs and resolution 1. These differ from manuscript settings and the later SOP branch. It cannot substitute for the missing upstream setup of the later full object. `min.features=2000` at object creation precedes the later >200-feature filter. Confirm metadata and the final parameter set before rerunning.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Cleaned code and execution

- [01_import_qc_sct.R](../scripts/snrna/01_import_qc_sct.R)
- [02_integrate.R](../scripts/snrna/02_integrate.R)

The cleaned import/QC/SCT script reconstructs the commented March upstream branch with required metadata/settings. The integration script implements its active tail. It can also start from the actual prepared Hm_merged object. This does not establish which branch produced the final paper object.

See [run commands and input formats](../scripts/README.md), [cleanup notes](../docs/CLEANUP_NOTES.md) and the [original source archive](../archive/README.md).

[Return to analysis index](../README.md#analysis-index).
