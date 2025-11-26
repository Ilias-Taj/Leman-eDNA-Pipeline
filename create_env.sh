#!/usr/bin/env bash
set -euo pipefail

# create_env.sh
# Create a conda/mamba environment in the './env' folder.
# Now uses the environment.yml file to install all tools.

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
else
  if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: neither mamba nor conda found on PATH. Please install Miniconda/Anaconda first." >&2
    exit 2
  fi
  echo "Using conda. Creating environment at: $ENV_PREFIX"
  conda env create -f "$YAML_FILE" -p "$ENV_PREFIX"
fi

echo ""
echo "Environment successfully created at: $ENV_PREFIX"
echo "Tools installed: cutadapt, filtlong, minimap2, kraken2, dorado"
echo "Activate it with: conda activate $ENV_PREFIX"