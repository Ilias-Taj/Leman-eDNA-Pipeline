#!/usr/bin/env bash
set -eo pipefail
# scripts/run_both_datasets.sh
# Runs the full eDNA pipeline sequentially: Water first, then Soil.
# Mirrors the argument style of run_full_pipeline.sh.
#
# Usage (foreground):
#   bash scripts/run_both_datasets.sh --db_18S pr2 --db_COI porter --db_JEDI pr2
#
# Usage (background):
#   nohup bash scripts/run_both_datasets.sh --db_18S pr2 --db_COI porter --db_JEDI pr2 \
#     > out/logs/pipeline_combined.log 2>&1 &
#   tail -f out/logs/pipeline_water.log

# Changelog:
#   2026-05-12  Refactor to accept CLI arguments (--db_18S, --db_COI, --db_JEDI,
#               --threads, --water_root, --soil_root) instead of env vars.
#               Add usage() help function.
# ── Defaults ──────────────────────────────────────────────────
DB_18S="pr2"
DB_COI="porter"
DB_JEDI="pr2"
THREADS=14
WATER_ROOT="data/Water_eDNA_18S_COI_14_01_26/fastq_pass"
SOIL_ROOT="data/Soil_eDNA_JEDI_COI_14_01_26/fastq_pass"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --db_18S  DB     18S database: pr2 (default), silva, or path to .udb
  --db_COI  DB     COI database: porter (default), midori2, ekoi, or path to .udb
  --db_JEDI DB     JEDI database: pr2 (default), silva, or path to .udb
  --threads N      Threads (default: $THREADS)
  --water_root DIR Water fastq_pass directory (default: $WATER_ROOT)
  --soil_root  DIR Soil fastq_pass directory  (default: $SOIL_ROOT)
  -h, --help       Show this help
EOF
}

# ── Argument parsing ──────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --db_18S)  DB_18S="$2";  shift 2;;
    --db_COI)  DB_COI="$2";  shift 2;;
    --db_JEDI) DB_JEDI="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --water_root) WATER_ROOT="$2"; shift 2;;
    --soil_root)  SOIL_ROOT="$2";  shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1"; usage; exit 1;;
  esac
done

mkdir -p out/logs
WATER_LOG="out/logs/pipeline_water.log"
SOIL_LOG="out/logs/pipeline_soil.log"

# ── Water dataset ─────────────────────────────────────────────
echo "=============================================="
echo "  Starting Water dataset (18S + COI)"
echo "  DB_18S=$DB_18S  DB_COI=$DB_COI"
echo "  $(date)"
echo "=============================================="

bash scripts/run_full_pipeline.sh \
  --root "$WATER_ROOT" \
  --markers 18S,COI \
  --threads "$THREADS" \
  --db_18S "$DB_18S" \
  --db_COI "$DB_COI" \
  2>&1 | tee "$WATER_LOG"

echo ""
echo "=============================================="
echo "  Water dataset finished at $(date)"
echo "=============================================="

# ── Soil dataset ──────────────────────────────────────────────
echo ""
echo "=============================================="
echo "  Starting Soil dataset (JEDI + COI)"
echo "  DB_JEDI=$DB_JEDI  DB_COI=$DB_COI"
echo "  $(date)"
echo "=============================================="

bash scripts/run_full_pipeline.sh \
  --root "$SOIL_ROOT" \
  --markers JEDI,COI \
  --threads "$THREADS" \
  --db_JEDI "$DB_JEDI" \
  --db_COI "$DB_COI" \
  2>&1 | tee "$SOIL_LOG"

echo ""
echo "=============================================="
echo "  Both datasets completed successfully!"
echo "  $(date)"
echo "=============================================="