#!/usr/bin/env bash
# Submit the real-data extraction/plot pipeline to Pegasus, never the login node.
set -euo pipefail
: "${TOOL:?Set absolute path to scsplice.py}"
: "${BAM:?Set BAM}"
: "${META:?Set metadata TSV}"
: "${SAMPLE:?Set sample identifier}"
: "${REGION:?Set contig:start-end}"
: "${PREFIX:?Set absolute output prefix}"
[[ $TOOL == /* && $BAM == /* && $META == /* && $PREFIX == /* ]] || {
  echo 'TOOL, BAM, META, PREFIX must be absolute paths' >&2; exit 2;
}
mkdir -p "$(dirname "$PREFIX")"
export TOOL BAM META SAMPLE REGION PREFIX
export STRAND="${STRAND:-forward}" ENV_NAME="${ENV_NAME:-scsplice}"
export STRIP_SUFFIX="${STRIP_SUFFIX:-0}" SPLICED_ONLY="${SPLICED_ONLY:-0}"
sbatch --parsable --partition="${PARTITION:-nano}" --cpus-per-task=1 \
  --mem="${MEMORY:-8G}" --time="${TIME:-00:30:00}" --job-name=scsplice \
  --output="${PREFIX}.%j.log" --export=ALL <<'SBATCH'
#!/usr/bin/env bash
set -euo pipefail
args=()
[[ $STRIP_SUFFIX != 1 ]] || args+=(--strip-barcode-suffix)
[[ $SPLICED_ONLY != 1 ]] || args+=(--spliced-only)
micromamba run -n "$ENV_NAME" python "$TOOL" extract \
  --bam "$BAM" --metadata "$META" --sample "$SAMPLE" --region "$REGION" \
  --strand "$STRAND" "${args[@]}" --output "$PREFIX"
micromamba run -n "$ENV_NAME" python "$TOOL" plot \
  --reads "${PREFIX}.reads.tsv" --region "$REGION" --output "${PREFIX}_plot"
SBATCH
