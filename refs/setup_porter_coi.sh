#!/usr/bin/env bash
set -euo pipefail
# Download and setup Porter CO1 Classifier v5.1.0 for VSEARCH SINTAX
# Source: https://github.com/terrimporter/CO1Classifier
# Contains 2.2M COI sequences from 236K taxa (185K species)
# Pre-built .udb ready for VSEARCH --sintax

cd "$(dirname "$0")/.."
echo "=== Downloading Porter CO1 v5.1.0 SINTAX database ==="
echo "Size: ~1.96 GB (compressed)"
echo ""

wget -c https://github.com/terrimporter/CO1Classifier/releases/download/SINTAX-COI-v5.1.0/SINTAX_COIv5.1.0.zip \
    -O refs/SINTAX_COIv5.1.0.zip

echo ""
echo "=== Extracting .udb ==="
unzip -o refs/SINTAX_COIv5.1.0.zip -d refs/
# The zip contains SINTAX_COIv5.1.0.udb
mv refs/SINTAX_COIv5.1.0.udb refs/porter_COI_v51.udb 2>/dev/null || true

echo ""
echo "=== Cleanup ==="
rm -f refs/SINTAX_COIv5.1.0.zip

echo ""
echo "Done! Database saved to: refs/porter_COI_v51.udb"
ls -lh refs/porter_COI_v51.udb
echo ""
echo "Usage:"
echo "  bash regenerate_taxonomy.sh --db_COI refs/porter_COI_v51.udb --dataset both"
