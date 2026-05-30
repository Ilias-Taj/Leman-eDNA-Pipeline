#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# ── PR2 version (update when new releases appear) ──
# Releases: https://github.com/pr2database/pr2database/releases
PR2_VERSION="${PR2_VERSION:-5.1.1}"

echo "=== Setting up PR2 v${PR2_VERSION} 18S database ==="

FASTA_GZ="refs/pr2_version_${PR2_VERSION}_SSU_dada2.fasta.gz"
SINTAX_FASTA="refs/pr2_18S_SINTAX.fasta"
UDB_FILE="refs/pr2_18S.udb"

# Download if not present
if [[ ! -f "$FASTA_GZ" ]]; then
    echo "Downloading PR2 v${PR2_VERSION} SSU DADA2 fasta (~52 MB)..."
    wget -c "https://github.com/pr2database/pr2database/releases/download/v${PR2_VERSION}/pr2_version_${PR2_VERSION}_SSU_dada2.fasta.gz" \
        -O "$FASTA_GZ"
else
    echo "PR2 fasta already present: $FASTA_GZ"
fi

echo ""
echo "=== Converting PR2 DADA2 format -> VSEARCH SINTAX format ==="
python3 refs/convert_pr2_to_sintax.py \
    --input "$FASTA_GZ" \
    --output "$SINTAX_FASTA"

echo ""
echo "=== Building VSEARCH .udb index ==="
./env/bin/vsearch --makeudb_usearch "$SINTAX_FASTA" \
    --output "$UDB_FILE"

echo ""
echo "=== Cleanup (keep .gz for reference) ==="
rm -f "$SINTAX_FASTA"

echo ""
echo "Done! Database saved to: $UDB_FILE"
ls -lh "$UDB_FILE"
echo ""
echo "Usage:"
echo "  bash scripts/regenerate_taxonomy.sh --db_18S $UDB_FILE --dataset water"
echo ""
echo "To use a different PR2 version:"
echo "  PR2_VERSION=5.2.0 bash refs/setup_pr2_18s.sh"