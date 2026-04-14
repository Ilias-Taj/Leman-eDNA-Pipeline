#!/usr/bin/env bash
set -eo pipefail
# regenerate_taxonomy.sh
# Re-runs taxonomy assignment (step 5) and comprehensive summary (step 7)
# for both Water and Soil datasets using the fixed confidence-filtered scripts.
#
# Usage:
#   bash regenerate_taxonomy.sh              # auto-detect (prefers eKOI)
#   bash regenerate_taxonomy.sh --db eKOI    # force eKOI
#   bash regenerate_taxonomy.sh --db midori2 # force MIDORI2

ENV_PREFIX="./env"
THREADS=14

export PATH="$ENV_PREFIX/bin:$PATH"

# Parse --db argument
DB_CHOICE=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --db)
      DB_CHOICE="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: bash regenerate_taxonomy.sh [--db eKOI|midori2]" >&2
      exit 1
      ;;
  esac
done

echo "=============================================="
echo "  Regenerating taxonomy (confidence-filtered)"
echo "  $(date)"
echo "=============================================="

# Determine COI database
if [ -n "$DB_CHOICE" ]; then
  case "$DB_CHOICE" in
    eKOI|ekoi)
      COI_DB="refs/eKOI_COI.udb"
      ;;
    midori2|MIDORI2|midori|MIDORI)
      COI_DB="refs/midori2_COI.udb"
      ;;
    *)
      echo "ERROR: Unknown database '$DB_CHOICE'. Use 'eKOI' or 'midori2'" >&2
      exit 1
      ;;
  esac
  if [ ! -f "$COI_DB" ]; then
    echo "ERROR: Database file not found: $COI_DB" >&2
    exit 1
  fi
else
  # Auto-detect: prefer eKOI, fall back to MIDORI2
  if [ -f refs/eKOI_COI.udb ]; then
    COI_DB="refs/eKOI_COI.udb"
  elif [ -f refs/midori2_COI.udb ]; then
    COI_DB="refs/midori2_COI.udb"
  else
    echo "ERROR: No COI database found in refs/" >&2
    exit 1
  fi
fi
echo "COI database: $COI_DB"

# -- Water dataset (18S + COI) --
WATER="out/Water_eDNA_18S_COI_14_01_26"
if [ -d "$WATER" ]; then
  echo ""
  echo "--- Water: taxonomy assignment ---"
  "$ENV_PREFIX/bin/python3" scripts/5_assign_taxonomy.py \
      --input_dir "$WATER" \
      --db_18S refs/silva_18s_v123.udb \
      --db_COI "$COI_DB" \
      --threads "$THREADS" \
      2>&1 | tee "$WATER/logs/taxonomy_assignment.log"

  echo ""
  echo "--- Water: comprehensive summary ---"
  "$ENV_PREFIX/bin/python3" scripts/7_comprehensive_taxonomy_summary.py \
      --input_dir "$WATER" \
      --db_18S refs/silva_18s_v123.udb \
      --db_COI "$COI_DB" \
      --skip_blast \
      2>&1 | tee "$WATER/logs/taxonomy_summary.log"
else
  echo "SKIP: $WATER not found"
fi

# -- Soil dataset (JEDI + COI) --
SOIL="out/Soil_eDNA_JEDI_COI_14_01_26"
if [ -d "$SOIL" ]; then
  echo ""
  echo "--- Soil: taxonomy assignment ---"
  "$ENV_PREFIX/bin/python3" scripts/5_assign_taxonomy.py \
      --input_dir "$SOIL" \
      --db_COI "$COI_DB" \
      --db_JEDI "$COI_DB" \
      --threads "$THREADS" \
      2>&1 | tee "$SOIL/logs/taxonomy_assignment.log"

  echo ""
  echo "--- Soil: comprehensive summary ---"
  "$ENV_PREFIX/bin/python3" scripts/7_comprehensive_taxonomy_summary.py \
      --input_dir "$SOIL" \
      --db_COI "$COI_DB" \
      --db_JEDI "$COI_DB" \
      --skip_blast \
      2>&1 | tee "$SOIL/logs/taxonomy_summary.log"
else
  echo "SKIP: $SOIL not found"
fi

echo ""
echo "=============================================="
echo "  Done! $(date)"
echo "=============================================="
