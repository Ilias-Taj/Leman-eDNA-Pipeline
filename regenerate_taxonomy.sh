#!/usr/bin/env bash
set -eo pipefail
# regenerate_taxonomy.sh
# Re-runs taxonomy assignment (step 5) and comprehensive summary (step 7)
# for both Water and Soil datasets using the fixed confidence-filtered scripts.
#
# Usage:
#   bash regenerate_taxonomy.sh                                    # auto-detect all
#   bash regenerate_taxonomy.sh --db_COI midori2 --db_JEDI eKOI   # per-marker
#   bash regenerate_taxonomy.sh --db_COI eKOI                     # only override COI
#
# Supported databases:
#   18S:  silva  (always refs/silva_18s_v123.udb)
#   COI:  eKOI | midori2
#   JEDI: eKOI | midori2

ENV_PREFIX="./env"
THREADS=14

export PATH="$ENV_PREFIX/bin:$PATH"

# Defaults (empty = auto-detect)
DB_COI_CHOICE=""
DB_JEDI_CHOICE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --db_COI)  DB_COI_CHOICE="$2";  shift 2 ;;
    --db_JEDI) DB_JEDI_CHOICE="$2"; shift 2 ;;
    --db)
      # Shorthand: set both COI and JEDI to the same DB
      DB_COI_CHOICE="$2"
      DB_JEDI_CHOICE="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: bash regenerate_taxonomy.sh [--db eKOI|midori2] [--db_COI eKOI|midori2] [--db_JEDI eKOI|midori2]" >&2
      exit 1
      ;;
  esac
done

resolve_coi_db() {
  local choice="$1"
  local label="$2"  # COI or JEDI (for messages)

  if [ -n "$choice" ]; then
    case "$choice" in
      eKOI|ekoi)       echo "refs/eKOI_COI.udb" ;;
      midori2|MIDORI2|midori|MIDORI) echo "refs/midori2_COI.udb" ;;
      *)
        echo "ERROR: Unknown database '$choice' for $label. Use 'eKOI' or 'midori2'" >&2
        exit 1
        ;;
    esac
  else
    # Auto-detect: prefer eKOI
    if [ -f refs/eKOI_COI.udb ]; then
      echo "refs/eKOI_COI.udb"
    elif [ -f refs/midori2_COI.udb ]; then
      echo "refs/midori2_COI.udb"
    else
      echo "ERROR: No COI database found in refs/ for $label" >&2
      exit 1
    fi
  fi
}

COI_DB=$(resolve_coi_db "$DB_COI_CHOICE" "COI")
JEDI_DB=$(resolve_coi_db "$DB_JEDI_CHOICE" "JEDI")

# Verify files exist
for db in "$COI_DB" "$JEDI_DB"; do
  if [ ! -f "$db" ]; then
    echo "ERROR: Database file not found: $db" >&2
    exit 1
  fi
done

echo "=============================================="
echo "  Regenerating taxonomy (confidence-filtered)"
echo "  $(date)"
echo "=============================================="
echo "  18S  database: refs/silva_18s_v123.udb"
echo "  COI  database: $COI_DB"
echo "  JEDI database: $JEDI_DB"
echo "=============================================="

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
      --db_JEDI "$JEDI_DB" \
      --threads "$THREADS" \
      2>&1 | tee "$SOIL/logs/taxonomy_assignment.log"

  echo ""
  echo "--- Soil: comprehensive summary ---"
  "$ENV_PREFIX/bin/python3" scripts/7_comprehensive_taxonomy_summary.py \
      --input_dir "$SOIL" \
      --db_COI "$COI_DB" \
      --db_JEDI "$JEDI_DB" \
      --skip_blast \
      2>&1 | tee "$SOIL/logs/taxonomy_summary.log"
else
  echo "SKIP: $SOIL not found"
fi

echo ""
echo "=============================================="
echo "  Done! $(date)"
echo "=============================================="
