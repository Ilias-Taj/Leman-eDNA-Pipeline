#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

FASTA_GZ="refs/MIDORI2_UNIQ_NUC_SP_GB269_CO1_SINTAX.fasta.gz"
FASTA="refs/MIDORI2_UNIQ_NUC_SP_GB269_CO1_SINTAX.fasta"

echo "Building MIDORI2 COI UDB database..."

if [[ -f "$FASTA_GZ" && ! -f "$FASTA" ]]; then
    echo "Decompressing $FASTA_GZ..."
    gunzip -k "$FASTA_GZ"
fi

if [[ ! -f "$FASTA" ]]; then
    echo "ERROR: $FASTA not found."
    echo "Download from: https://reference-midori.info/download.php"
    echo "  File: MIDORI2_UNIQ_NUC_SP_GB269_CO1_SINTAX.fasta.gz"
    exit 1
fi

echo "This will take a few minutes for 3.1M sequences."
rm -f refs/midori2_COI.udb
./env/bin/vsearch --makeudb_usearch "$FASTA" --output refs/midori2_COI.udb
echo ""
echo "Done! Database saved to: refs/midori2_COI.udb"
ls -lh refs/midori2_COI.udb
