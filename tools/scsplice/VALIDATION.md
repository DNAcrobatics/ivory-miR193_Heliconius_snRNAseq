# Validation

Executed on 2026-10-21 with Python 3.12.14, pysam 0.24.1 and Matplotlib 3.10.8.

`python -m unittest -v test_scsplice.py`: 3 tests passed, including the complete
synthetic BAM → metadata matching → read table → two PDF plotting modes flow.

`bash -n submit_pegasus.sh`: passed.

The R metadata exporter and archived R scripts were not executed (R is unavailable
in this runtime). No biological analysis, alignment, clustering, real-data BAM
scan, cluster submission or published-figure reproduction was performed.

Source copies are unmodified. SOURCE_SHA256SUMS.txt records their checksums and
the shared R helper. The README records the GitHub source snapshot and paths.
