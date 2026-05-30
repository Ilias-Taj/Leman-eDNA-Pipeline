#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

echo "Setting up eKOI COI database..."

# Download raw FASTA from GitHub
if [[ ! -f refs/eKOI_taxonomy_PR2_ver1.fasta ]]; then
    echo "Downloading eKOI PR2 FASTA..."
    wget -O refs/eKOI_taxonomy_PR2_ver1.fasta \
        https://raw.githubusercontent.com/rubenmiguens/eKOI_taxonomy_database/main/eKOI_taxonomy_PR2_ver1.fasta
else
    echo "eKOI raw FASTA already present, skipping download."
fi

# Convert to SINTAX format
echo "Converting to SINTAX format..."
python3 refs/convert_ekoi_to_sintax.py \
    --input refs/eKOI_taxonomy_PR2_ver1.fasta \
    --output refs/eKOI_COI_SINTAX.fasta

# Build UDB
echo "Building UDB database..."
./env/bin/vsearch --makeudb_usearch refs/eKOI_COI_SINTAX.fasta --output refs/eKOI_COI.udb

echo ""
echo "Done! Database saved to: refs/eKOI_COI.udb"
ls -lh refs/eKOI_COI.udb
