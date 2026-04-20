#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
echo "Building MIDORI2 COI UDB database..."
echo "This will take a few minutes for 3.1M sequences."
rm -f refs/midori2_COI.udb
./env/bin/vsearch --makeudb_usearch refs/MIDORI2_UNIQ_NUC_SP_GB269_CO1_SINTAX.fasta --output refs/midori2_COI.udb
echo ""
echo "Done! Database saved to: refs/midori2_COI.udb"
ls -lh refs/midori2_COI.udb
