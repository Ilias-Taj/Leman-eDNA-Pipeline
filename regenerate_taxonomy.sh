#!/usr/bin/env bash
set -eo pipefail
# regenerate_taxonomy.sh
# Re-runs taxonomy assignment (step 5) and comprehensive summary (step 7)
# for both Water and Soil datasets using the fixed confidence-filtered scripts.
#
# Usage:
#   bash regenerate_taxonomy.sh                                              # both datasets, auto-detect
#   bash regenerate_taxonomy.sh --dataset water --db_18S silva --db_COI ekoi # water only, custom DBs
#   bash regenerate_taxonomy.sh --dataset soil --db ekoi                     # soil only, eKOI for COI+JEDI
#   bash regenerate_taxonomy.sh --db_18S silva --db_COI midori2 --db_JEDI eKOI
#   bash regenerate_taxonomy.sh --db midori2                                 # set COI+JEDI
#
# Supported databases:
#   --db_18S:  silva (default: refs/silva_18s_v123.udb)
#              Add new 18S .udb files to refs/ and extend resolve_18s_db()
#   --db_COI:  eKOI | midori2
#   --db_JEDI: eKOI | midori2 | pr2 (recommended for 18S-like JEDI marker)
#   --db:      shorthand to set COI + JEDI to the same database

ENV_PREFIX="./env"
THREADS=14

export PATH="$ENV_PREFIX/bin:$PATH"

# Defaults (empty = auto-detect)
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
    --db)
      # Shorthand: set COI and JEDI to the same DB
      DB_COI_CHOICE="$2"
      DB_JEDI_CHOICE="$2"
      shift 2
      ;;
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
      echo "Usage: bash regenerate_taxonomy.sh [--dataset water|soil|both] [--db_18S DB] [--db_COI DB] [--db_JEDI DB] [--db DB]" >&2
      exit 1
      ;;
  esac
done

resolve_18s_db() {
  local choice="$1"

  if [ -n "$choice" ]; then
    case "$choice" in
      silva|SILVA)
        echo "refs/silva_18s_v123.udb" ;;
      pr2|PR2) echo "refs/pr2_18S_v511.udb" ;;
      *)
        # Treat as direct path if it exists
        if [ -f "$choice" ]; then
          echo "$choice"
        else
          echo "ERROR: Unknown 18S database '$choice'. Use 'silva' or a path to a .udb file" >&2
          exit 1
        fi
        ;;
    esac
  else
    # Default to SILVA
    if [ -f refs/silva_18s_v123.udb ]; then
      echo "refs/silva_18s_v123.udb"
    else
      echo "ERROR: No 18S database found (refs/silva_18s_v123.udb)" >&2
      exit 1
    fi
  fi
}

resolve_coi_db() {
  local choice="$1"
  local label="$2"  # COI or JEDI (for messages)

  if [ -n "$choice" ]; then
    case "$choice" in
      eKOI|ekoi)       echo "refs/eKOI_COI.udb" ;;
      midori2|MIDORI2|midori|MIDORI) echo "refs/midori2_COI.udb" ;;
      porter|PORTER|porter_coi) echo "refs/porter_COI_v51.udb" ;;
      pr2|PR2) echo "refs/pr2_18S_v511.udb" ;;
      *)
        # Treat as direct path if it exists
        if [ -f "$choice" ]; then
          echo "$choice"
        else
          echo "ERROR: Unknown database '$choice' for $label. Use 'eKOI', 'midori2', 'porter', 'pr2', or a path to a .udb file" >&2
          exit 1
        fi
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

DB_18S=$(resolve_18s_db "$DB_18S_CHOICE")
COI_DB=$(resolve_coi_db "$DB_COI_CHOICE" "COI")
JEDI_DB=$(resolve_coi_db "$DB_JEDI_CHOICE" "JEDI")

# Verify files exist
for db in "$DB_18S" "$COI_DB" "$JEDI_DB"; do
  if [ ! -f "$db" ]; then
    echo "ERROR: Database file not found: $db" >&2
    exit 1
  fi
done

echo "=============================================="
echo "  Regenerating taxonomy (confidence-filtered)"
echo "  $(date)"
echo "=============================================="
echo "  18S  database: $DB_18S"
echo "  COI  database: $COI_DB"
echo "  JEDI database: $JEDI_DB"
echo "=============================================="
echo "  Dataset(s):    $DATASET"
echo "  (outputs will be placed under taxonomy/<tag>/ and taxonomy_summary/<tag>/"
echo "   where <tag> is auto-derived from the DB filenames by scripts 5 and 7)"
echo "=============================================="

# -- Water dataset (18S + COI) --
WATER="out/Water_eDNA_18S_COI_14_01_26"
if [[ "$DATASET" == "water" || "$DATASET" == "both" ]] && [ -d "$WATER" ]; then
  echo ""
  echo "--- Water: taxonomy assignment ---"
  "$ENV_PREFIX/bin/python3" scripts/5_assign_taxonomy.py \
      --input_dir "$WATER" \
      --db_18S "$DB_18S" \
      --db_COI "$COI_DB" \
      --threads "$THREADS" \
      2>&1 | tee "$WATER/logs/taxonomy_assignment.log"

  echo ""
  echo "--- Water: comprehensive summary ---"
  "$ENV_PREFIX/bin/python3" scripts/7_comprehensive_taxonomy_summary.py \
      --input_dir "$WATER" \
      --db_18S "$DB_18S" \
      --db_COI "$COI_DB" \
      --skip_blast \
      2>&1 | tee "$WATER/logs/taxonomy_summary.log"
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
