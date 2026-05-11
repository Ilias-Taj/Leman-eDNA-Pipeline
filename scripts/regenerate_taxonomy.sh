#!/usr/bin/env bash
set -eo pipefail

# Ensure we run from the repo root regardless of where the script is called from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

# scripts/regenerate_taxonomy.sh
# Re-runs taxonomy assignment (step 5) and comprehensive summary (step 7).
# Only processes markers for which a database is explicitly provided.
# Markers without a --db_* argument are skipped entirely.
#
# Usage:
#   bash scripts/regenerate_taxonomy.sh --db_18S silva --db_COI midori2
#   bash scripts/regenerate_taxonomy.sh --dataset water --db_COI ekoi
#   bash scripts/regenerate_taxonomy.sh --db_18S pr2 --db_COI porter --db_JEDI pr2
#   bash scripts/regenerate_taxonomy.sh --db midori2          # shorthand for --db_COI
#
# Supported databases:
#   --db_18S:  silva | pr2
#   --db_COI:  midori2 | ekoi | porter
#   --db_JEDI: silva | pr2  (JEDI targets rRNA V4-V5, same as 18S)
#   --db:      shorthand for --db_COI

ENV_PREFIX="./env"
THREADS=14

export PATH="$ENV_PREFIX/bin:$PATH"

# Defaults: empty = skip marker
DB_18S_CHOICE=""
DB_COI_CHOICE=""
DB_JEDI_CHOICE=""
DATASET="both"  # water | soil | both

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --db_18S)  DB_18S_CHOICE="$2";  shift 2 ;;
    --db_COI)  DB_COI_CHOICE="$2";  shift 2 ;;
    --db_JEDI) DB_JEDI_CHOICE="$2"; shift 2 ;;
    --db)      DB_COI_CHOICE="$2";  shift 2 ;;
    --dataset)
      DATASET="$2"
      case "$DATASET" in
        water|soil|both) ;;
        *) echo "ERROR: --dataset must be water|soil|both (got: $DATASET)" >&2; exit 1 ;;
      esac
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      echo "Usage: bash scripts/regenerate_taxonomy.sh [--dataset water|soil|both] [--db_18S DB] [--db_COI DB] [--db_JEDI DB] [--db DB]" >&2
      exit 1
      ;;
  esac
done

# Resolve a DB name to a .udb path. Returns empty string if choice is empty.
resolve_db() {
  local choice="$1"
  local label="$2"  # for error messages

  [ -z "$choice" ] && return 0

  case "$choice" in
    silva|SILVA)                      echo "refs/silva_18s_v123.udb" ;;
    pr2|PR2)                          echo "refs/pr2_18S_v511.udb" ;;
    midori2|MIDORI2|midori|MIDORI)    echo "refs/midori2_COI.udb" ;;
    eKOI|ekoi)                        echo "refs/eKOI_COI.udb" ;;
    porter|PORTER|porter_coi)         echo "refs/porter_COI_v51.udb" ;;
    *)
      if [ -f "$choice" ]; then
        echo "$choice"
      else
        echo "ERROR: Unknown database '$choice' for $label. Use silva|pr2|midori2|ekoi|porter or a .udb path" >&2
        exit 1
      fi
      ;;
  esac
}

DB_18S=$(resolve_db "$DB_18S_CHOICE" "18S")
COI_DB=$(resolve_db "$DB_COI_CHOICE" "COI")
JEDI_DB=$(resolve_db "$DB_JEDI_CHOICE" "JEDI")

# Verify chosen DBs exist
for db in "$DB_18S" "$COI_DB" "$JEDI_DB"; do
  if [ -n "$db" ] && [ ! -f "$db" ]; then
    echo "ERROR: Database file not found: $db" >&2
    exit 1
  fi
done

# Check at least one DB was provided
if [ -z "$DB_18S" ] && [ -z "$COI_DB" ] && [ -z "$JEDI_DB" ]; then
  echo "ERROR: No database specified. Provide at least one of --db_18S, --db_COI, --db_JEDI, or --db" >&2
  exit 1
fi

echo "=============================================="
echo "  Regenerating taxonomy"
echo "  $(date)"
echo "=============================================="
[ -n "$DB_18S" ]  && echo "  18S  database: $DB_18S"  || echo "  18S  database: (skip)"
[ -n "$COI_DB" ]  && echo "  COI  database: $COI_DB"  || echo "  COI  database: (skip)"
[ -n "$JEDI_DB" ] && echo "  JEDI database: $JEDI_DB" || echo "  JEDI database: (skip)"
echo "  Dataset(s):    $DATASET"
echo "  Output layout: taxonomy/{MARKER}/{db}/"
echo "=============================================="

# Build --db_* args only for chosen markers
build_db_args() {
  local dataset="$1"  # water or soil
  local args=""
  if [ "$dataset" = "water" ]; then
    [ -n "$DB_18S" ]  && args="$args --db_18S $DB_18S"
    [ -n "$COI_DB" ]  && args="$args --db_COI $COI_DB"
  elif [ "$dataset" = "soil" ]; then
    [ -n "$COI_DB" ]  && args="$args --db_COI $COI_DB"
    [ -n "$JEDI_DB" ] && args="$args --db_JEDI $JEDI_DB"
  fi
  echo "$args"
}

# -- Water dataset (18S + COI) --
WATER="out/Water_eDNA_18S_COI_14_01_26"
if [[ "$DATASET" == "water" || "$DATASET" == "both" ]] && [ -d "$WATER" ]; then
  WATER_ARGS=$(build_db_args water)
  if [ -n "$WATER_ARGS" ]; then
    echo ""
    echo "--- Water: taxonomy assignment ---"
    "$ENV_PREFIX/bin/python3" scripts/5_assign_taxonomy.py \
        --input_dir "$WATER" \
        $WATER_ARGS \
        --threads "$THREADS" \
        2>&1 | tee "$WATER/logs/taxonomy_assignment.log"

    echo ""
    echo "--- Water: comprehensive summary ---"
    "$ENV_PREFIX/bin/python3" scripts/7_comprehensive_taxonomy_summary.py \
        --input_dir "$WATER" \
        $WATER_ARGS \
        --skip_blast \
        2>&1 | tee "$WATER/logs/taxonomy_summary.log"
  else
    echo ""
    echo "[SKIP] Water: no relevant DB provided (need --db_18S or --db_COI)"
  fi
else
  if [[ "$DATASET" == "soil" ]]; then
    echo "[dataset=soil] Skipping Water"
  else
    echo "SKIP: $WATER not found"
  fi
fi

# -- Soil dataset (JEDI + COI) --
SOIL="out/Soil_eDNA_JEDI_COI_14_01_26"
if [[ "$DATASET" == "soil" || "$DATASET" == "both" ]] && [ -d "$SOIL" ]; then
  SOIL_ARGS=$(build_db_args soil)
  if [ -n "$SOIL_ARGS" ]; then
    echo ""
    echo "--- Soil: taxonomy assignment ---"
    "$ENV_PREFIX/bin/python3" scripts/5_assign_taxonomy.py \
        --input_dir "$SOIL" \
        $SOIL_ARGS \
        --threads "$THREADS" \
        2>&1 | tee "$SOIL/logs/taxonomy_assignment.log"

    echo ""
    echo "--- Soil: comprehensive summary ---"
    "$ENV_PREFIX/bin/python3" scripts/7_comprehensive_taxonomy_summary.py \
        --input_dir "$SOIL" \
        $SOIL_ARGS \
        --skip_blast \
        2>&1 | tee "$SOIL/logs/taxonomy_summary.log"
  else
    echo ""
    echo "[SKIP] Soil: no relevant DB provided (need --db_COI or --db_JEDI)"
  fi
else
  if [[ "$DATASET" == "water" ]]; then
    echo "[dataset=water] Skipping Soil"
  else
    echo "SKIP: $SOIL not found"
  fi
fi

echo ""
echo "=============================================="
echo "  Done! $(date)"
echo "=============================================="
