# miR-193 result validation

The source is intake commit `4adaff5b3b2143542a409df2d63261f7718bfa4d` on `mir-intake-20261008_203210`. Its 14 archive entries are retained as the exact Git blobs supplied by the user.

## Checks executed

- All 12 available text files match their entries in `SHA256SUMS`.
- Parsed 2,079 verbose predictions; grouped counts and minimum energies agree exactly with all 1,265 rows of the supplied gene summary.
- The reconstructed predicted-gene set equals the supplied list.
- The comparative join agrees with all 733 Table S8 hit counts and Boolean annotations and recovers 15 conserved candidate rows.
- The only energy discrepancy is KrT95D / Hmel.Hmel208001oG99: missing in S8, −26.4 in archived predictions. The original S8 is retained unchanged.
- Empty, activation-error and truncated prediction input are rejected by the parser.
- Bash syntax passes. Mock launcher execution verifies micromamba argument order, retained prediction flags, overwrite rejection and failure handling. Mock execution does not validate RNAhybrid predictions.

The compressed UTR FASTA is retained through its original Git blob (`c877c872e72f2ee9f6bb549ff9837a3646aa84fd`). Its SHA-256 checksum is supplied in the intake manifest but has not been independently checked here because the connector could not retrieve the binary. Neither its FASTA content nor the longest-transcript selection was inspected. The launcher accepts gzip input, but gzip execution was not included in the mock checks.

## Provenance limits

RNAhybrid was not rerun. All three supplied job logs contain shell-activation errors. The cleaned launcher uses `micromamba run`, retaining the supplied `-s 3utr_fly -m 10000 -b 10 -e -25`. Manuscript version 2.1.2 and cutoff −20 remain reported values, not verified execution settings.

Gene collapsing and comparative joins are reconstructed procedures validated against the supplied outputs, not recovered historical parser code. The conserved-candidate rule is predicted hit plus DE annotation in at least two studies. Primary cross-species inputs and orthology/filtering provenance are still required for independent reproduction of those flags.
