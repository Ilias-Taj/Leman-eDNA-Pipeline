#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# ── eKOI version (update when new versions appear) ──
# Source: https://github.com/rubenmiguens/eKOI_taxonomy_database
EKOI_VERSION="${EKOI_VERSION:-1}"

echo "=== Setting up eKOI COI database (ver${EKOI_VERSION}) ==="

RAW_FASTA="refs/eKOI_taxonomy_PR2_ver${EKOI_VERSION}.fasta"
SINTAX_FASTA="refs/eKOI_COI_SINTAX.fasta"
UDB_FILE="refs/eKOI_COI.udb"

# Download raw FASTA from GitHub
if [[ ! -f "$RAW_FASTA" ]]; then
    echo "Downloading eKOI PR2 ver${EKOI_VERSION} FASTA..."
    wget -O "$RAW_FASTA" \
        "https://raw.githubusercontent.com/rubenmiguens/eKOI_taxonomy_database/main/eKOI_taxonomy_PR2_ver${EKOI_VERSION}.fasta"
else
    echo "eKOI raw FASTA already present: $RAW_FASTA"
fi

# Convert to SINTAX format
echo "Converting to SINTAX format..."
python3 refs/convert_ekoi_to_sintax.py \
    --input "$RAW_FASTA" \
    --output "$SINTAX_FASTA"

# Build UDB
echo "Building UDB database..."
./env/bin/vsearch --makeudb_usearch "$SINTAX_FASTA" --output "$UDB_FILE"

echo ""
echo "Done! Database saved to: $UDB_FILE"
ls -lh "$UDB_FILE"
echo ""
echo "Usage:"
echo "  bash scripts/regenerate_taxonomy.sh --db_COI $UDB_FILE --dataset both"
echo ""
echo "To use a different eKOI version:"
echo "  EKOI_VERSION=2 bash refs/setup_ekoi_coi.sh"