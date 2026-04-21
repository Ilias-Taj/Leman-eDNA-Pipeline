#!/usr/bin/env bash
set -eo pipefail
# scripts/run_both_datasets.sh
# Runs the full eDNA pipeline sequentially: Water first, then Soil.
# Each dataset gets its own independent log file while output is also
# printed to the terminal so you can follow progress live.
#
# Usage (foreground – see output live):
#   bash scripts/run_both_datasets.sh
#
# Usage (background – follow with tail):
#   nohup bash scripts/run_both_datasets.sh > out/logs/pipeline_combined.log 2>&1 &
#   tail -f out/logs/pipeline_water.log   # or pipeline_soil.log

WATER_LOG="out/logs/pipeline_water.log"
SOIL_LOG="out/logs/pipeline_soil.log"
mkdir -p out/logs

# ── Water dataset ─────────────────────────────────────────────
echo "=============================================="
echo "  Starting Water dataset (18S + COI)"
echo "  $(date)"
echo "=============================================="

bash scripts/run_full_pipeline.sh \
  --root data/Water_eDNA_18S_COI_14_01_26/fastq_pass \
  --markers 18S,COI \
  --threads 14 2>&1 | tee "$WATER_LOG"

echo ""
echo "=============================================="
echo "  Water dataset finished at $(date)"
echo "=============================================="

# ── Soil dataset ──────────────────────────────────────────────
echo ""
echo "=============================================="
echo "  Starting Soil dataset (JEDI + COI)"
echo "  $(date)"
echo "=============================================="

bash scripts/run_full_pipeline.sh \
  --root data/Soil_eDNA_JEDI_COI_14_01_26/fastq_pass \
  --markers JEDI,COI \
  --threads 14 2>&1 | tee "$SOIL_LOG"

echo ""
echo "=============================================="
echo "  Both datasets completed successfully!"
echo "  $(date)"
echo "=============================================="
