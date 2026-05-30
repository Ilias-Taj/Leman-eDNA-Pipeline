#!/usr/bin/env bash
set -euo pipefail
# Download and setup Porter CO1 Classifier for VSEARCH SINTAX
# Source: https://github.com/terrimporter/CO1Classifier
# Contains 2.2M COI sequences from 236K taxa (185K species)

cd "$(dirname "$0")/.."

# ── Porter COI version (update when new releases appear) ──
# Releases: https://github.com/terrimporter/CO1Classifier/releases
PORTER_VERSION="${PORTER_VERSION:-5.1.0}"
PORTER_SHORT="${PORTER_VERSION%.*}${PORTER_VERSION##*.}"  # 5.1.0 -> 510 (for filenames)

echo "=== Downloading Porter CO1 v${PORTER_VERSION} SINTAX database ==="
echo "Size: ~1.96 GB (compressed)"
echo ""

ZIP_FILE="refs/SINTAX_COIv${PORTER_VERSION}.zip"
UDB_FILE="refs/porter_COI_v${PORTER_SHORT%0}.udb"

if [[ ! -f "$ZIP_FILE" ]]; then
    wget -c "https://github.com/terrimporter/CO1Classifier/releases/download/SINTAX-COI-v${PORTER_VERSION}/SINTAX_COIv${PORTER_VERSION}.zip" \
        -O "$ZIP_FILE"
else
    echo "ZIP already present: $ZIP_FILE"
fi

echo ""
echo "=== Extracting .udb ==="
unzip -o "$ZIP_FILE" -d refs/
# The zip contains trained/sintax.udb inside
mv refs/trained/sintax.udb "$UDB_FILE"
rm -rf refs/trained

echo ""
echo "=== Cleanup ==="
rm -f "$ZIP_FILE"

echo ""
echo "Done! Database saved to: $UDB_FILE"
ls -lh "$UDB_FILE"
echo ""
echo "Usage:"
echo "  bash scripts/regenerate_taxonomy.sh --db_COI $UDB_FILE --dataset both"
echo ""
echo "To use a different Porter version:"
echo "  PORTER_VERSION=5.2.0 bash refs/setup_porter_coi.sh"