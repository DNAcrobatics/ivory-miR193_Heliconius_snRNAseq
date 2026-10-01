# Candidate miR-193 target prediction

**Paper:** Figure 7 target highlighting.

**Status:** Gene list supplied; prediction code needed.

Recover the predicted target set and its joins to expression comparisons.

## Inputs

- Mature miR-193 FASTA, longest 3′ UTR FASTA and their derivation scripts.
- RNAhybrid command/job files and raw prediction outputs.
- Complete `with_miRhit` gene set and final target table.

## Outputs

- Reproducible candidate predictions and annotated DE tables for target highlighting.

## Software

RNAhybrid 2.1.2 is stated in the methods, with 3utr_fly and energy cutoff -20; executed commands and version logs still needed. R for highlighting.

## Execution and interpretation

The supplied fragment defines a selected 15-gene list, not the prediction pipeline or complete hit set. Recover the longest-UTR selection rules, identifier mapping, mature sequence and filtering. Preserve the distinction between all predicted hits and this selected subset.

See [intake review](../docs/INTAKE_REVIEW.md) for unresolved decisions and [file requests](../docs/FILE_REQUESTS.md) for the complete acquisition checklist.

## Code

### Recovered selected target list

**Source:** [ivory_reanalysis.expanded.pasted.r, lines 1434–1449](../archive/original/ivory_proj_transfer_aug2025/ivory_reanalysis.expanded.pasted.r#L1434-L1449).

**Archive treatment:** exact excerpt; no syntax repairs, parameter changes, or inferred session setup. This is a recovered session fragment, not a standalone script.

```r
## genes that have a mir hit and are expressed in other species in notable ways
top15_mir_multi <- c("Hmel.Hmel220003oG243",
"Hmel.Hmel221001oG385",
"Hmel.Hmel217001oG324",
"Hmel.Hmel201001oG305",
"Hmel.Hmel204001oG141",
"Hmel.Hmel220003oG261",
"Hmel.Hmel206001oG635",
"Hmel.Hmel213001oG325",
"Hmel.Hmel205001oG141",
"Hmel.Hmel214004oG620",
"Hmel.Hmel213001oG409",
"Hmel.Hmel212001oG92",
"Hmel.Hmel201001oG657",
"Hmel.Hmel207001oG198",
"Hmel.Hmel203003oG386")
```


[Return to analysis index](../README.md#analysis-index).
