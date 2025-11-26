#!/usr/bin/env bash
# Build a concatenated SANGER reference FASTA from the dataset
# Usage: ./scripts/1_build_sanger_db.sh

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$ROOT_DIR/raw_reads"
OUT_DIR="$ROOT_DIR/refs"
OUT_FA="$OUT_DIR/sanger_db.fa"

mkdir -p "$OUT_DIR"

echo "Building SANGER reference DB at: $OUT_FA"

>"$OUT_FA"
find "$DATA_DIR" -type f -path "*/SANGER/*.fasta" | while read -r f; do
  species_dir=$(basename "$(dirname "$(dirname "$f")")")
  awk -v prefix="$species_dir" 'BEGIN{RS=">"; ORS=""} NR>1{split($0,lines,"\n"); hdr=lines[1]; seq=""; for(i=2;i<=length(lines);i++){seq=seq lines[i]"\n"} print ">"prefix"|"hdr"\n"seq}' "$f" >> "$OUT_FA"
done

echo "Done. Created $OUT_FA"
