#!/usr/bin/env bash
set -euo pipefail
# Download PR2 v5.1.1 SSU database and convert to VSEARCH SINTAX format
# Source: https://github.com/pr2database/pr2database
# Contains ~240K 18S rRNA sequences with expert-curated taxonomy
# Excellent for protists, ciliates, diatoms, plankton

cd "$(dirname "$0")/.."
echo "=== Downloading PR2 v5.1.1 SSU DADA2 fasta ==="
echo "Size: ~52 MB (compressed)"
echo ""

wget -c https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_dada2.fasta.gz \
    -O refs/pr2_version_5.1.1_SSU_dada2.fasta.gz

echo ""
echo "=== Converting PR2 DADA2 format -> VSEARCH SINTAX format ==="
python refs/convert_pr2_to_sintax.py \
    --input refs/pr2_version_5.1.1_SSU_dada2.fasta.gz \
    --output refs/pr2_18S_SINTAX.fasta

echo ""
echo "=== Building VSEARCH .udb index ==="
./env/bin/vsearch --makeudb_usearch refs/pr2_18S_SINTAX.fasta \
    --output refs/pr2_18S_v511.udb

echo ""
echo "=== Cleanup (keep .gz for reference) ==="
rm -f refs/pr2_18S_SINTAX.fasta

echo ""
echo "Done! Database saved to: refs/pr2_18S_v511.udb"
ls -lh refs/pr2_18S_v511.udb
echo ""
echo "Usage:"
echo "  bash regenerate_taxonomy.sh --db_18S refs/pr2_18S_v511.udb --dataset water"
