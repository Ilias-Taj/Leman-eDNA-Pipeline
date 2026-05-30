#!/usr/bin/env bash
set -euo pipefail

# create_env.sh — Create conda/mamba environment for the eDNA pipeline
# Installs all dependencies defined in environment.yml

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_PREFIX="$ROOT_DIR/env"
YAML_FILE="$ROOT_DIR/environment.yml"

if [ ! -f "$YAML_FILE" ]; then
    echo "ERROR: environment.yml not found at: $YAML_FILE" >&2
    exit 1
fi

if command -v mamba >/dev/null 2>&1; then
  echo "Found mamba. Creating environment at: $ENV_PREFIX"
  mamba env create -f "$YAML_FILE" -p "$ENV_PREFIX"
elif command -v conda >/dev/null 2>&1; then
  echo "Using conda. Creating environment at: $ENV_PREFIX"
  conda env create -f "$YAML_FILE" -p "$ENV_PREFIX"
else
  echo "ERROR: neither mamba nor conda found on PATH." >&2
  echo "Install Miniconda: https://docs.conda.io/en/latest/miniconda.html" >&2
  exit 2
fi

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo " Environment created at: $ENV_PREFIX"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo "Conda packages installed:"
echo "  filtlong, vsearch, samtools, cutadapt, biopython, pandas,"
echo "  matplotlib, seaborn, numpy, ipykernel, jupyter"
echo ""
echo "Pip packages installed:"
echo "  spoa (SPOA consensus generation)"
echo ""
echo "Activate with:"
echo "  conda activate $ENV_PREFIX"
echo ""
echo "═══════════════════════════════════════════════════════════════"
echo " External tools (install manually in tools/):"
echo "═══════════════════════════════════════════════════════════════"
echo "  minimap2 v2.30+  — https://github.com/lh3/minimap2/releases"
echo "                     Download pre-compiled binary for your platform"
echo "                     (newer versions are backward-compatible)"
echo ""
echo "  isONclust3       — https://github.com/aljpetri/isONclust3"
echo "                     cargo build --release (requires Rust toolchain)"
echo "                     or download a release binary if available"
echo ""
echo "  Place binaries in tools/ directory. Newer versions generally work."
echo "═══════════════════════════════════════════════════════════════"
