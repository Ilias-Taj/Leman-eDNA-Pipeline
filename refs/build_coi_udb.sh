#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# ── MIDORI2 COI version (update when new GenBank releases appear) ──
# Download from: https://reference-midori.info/download.php
# Look for UNIQ_NUC_SP file with SINTAX format
MIDORI_GB_VERSION="${MIDORI_GB_VERSION:-}"

echo "Building MIDORI2 COI UDB database..."

# Auto-detect: find any MIDORI2 SP SINTAX fasta in refs/
if [[ -z "$MIDORI_GB_VERSION" ]]; then
    FASTA=$(find refs/ -maxdepth 1 -name "MIDORI2_UNIQ_NUC_SP_GB*_CO1_SINTAX.fasta" 2>/dev/null | sort -V | tail -1 || true)
    FASTA_GZ=$(find refs/ -maxdepth 1 -name "MIDORI2_UNIQ_NUC_SP_GB*_CO1_SINTAX.fasta.gz" 2>/dev/null | sort -V | tail -1 || true)
else
    FASTA="refs/MIDORI2_UNIQ_NUC_SP_GB${MIDORI_GB_VERSION}_CO1_SINTAX.fasta"
    FASTA_GZ="refs/MIDORI2_UNIQ_NUC_SP_GB${MIDORI_GB_VERSION}_CO1_SINTAX.fasta.gz"
fi

# Decompress if needed
if [[ -n "$FASTA_GZ" && -f "$FASTA_GZ" && -z "$FASTA" ]]; then
    echo "Decompressing $FASTA_GZ..."
    gunzip -k "$FASTA_GZ"
    FASTA="${FASTA_GZ%.gz}"
fi

if [[ -z "$FASTA" || ! -f "$FASTA" ]]; then
    echo "ERROR: No MIDORI2 SINTAX fasta found in refs/"
    echo "Download from: https://reference-midori.info/download.php"
    echo "  File: MIDORI2_UNIQ_NUC_SP_GB<version>_CO1_SINTAX.fasta.gz"
    echo "  Place in: refs/"
    echo ""
    echo "Or set version explicitly: MIDORI_GB_VERSION=269 bash refs/build_coi_udb.sh"
    exit 1
fi

echo "Using: $FASTA"
echo "This will take a few minutes for 3.1M sequences."
rm -f refs/midori2_COI.udb
./env/bin/vsearch --makeudb_usearch "$FASTA" --output refs/midori2_COI.udb
echo ""
echo "Done! Database saved to: refs/midori2_COI.udb"
ls -lh refs/midori2_COI.udb