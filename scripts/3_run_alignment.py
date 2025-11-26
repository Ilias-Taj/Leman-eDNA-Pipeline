#!/usr/bin/env python3
"""
3_run_alignment.py

Align preprocessed reads against a targeted reference database using Minimap2.

Usage example:
    python 3_run_alignment.py --reads_file out/filtered.fastq --db_file test_db.fa --output_dir out --threads 4

This script will produce:
  - <output_dir>/alignments.paf
"""

import argparse
import subprocess
import sys
from pathlib import Path

def run_cmd(cmd):
    """Prints and runs a command, exiting on failure."""
    print("Running:", " ".join(map(str, cmd)))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)
    except FileNotFoundError as e:
        print(f"ERROR: Command not found. Is {cmd[0]} installed and in your PATH?", file=sys.stderr)
        print(e, file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Alignment: Minimap2")
    parser.add_argument("--reads_file", required=True, help="Path to filtered FASTQ file (from 2_run_preprocessing.py)")
    parser.add_argument("--db_file", required=True, help="Path to reference FASTA database (e.g., MiFish_12S.fa)")
    parser.add_argument("--output_dir", required=True, help="Directory to save alignment file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")

    args = parser.parse_args()

    reads_path = Path(args.reads_file)
    db_path = Path(args.db_file)
    out_dir = Path(args.output_dir)
    threads = args.threads

    if not reads_path.exists():
        print(f"ERROR: input reads file does not exist: {reads_path}", file=sys.stderr)
        sys.exit(2)
    if not db_path.exists():
        print(f"ERROR: database file does not exist: {db_path}", file=sys.stderr)
        sys.exit(2)

    out_dir.mkdir(parents=True, exist_ok=True)
    paf_output_path = out_dir / "alignments.paf"

    minimap_cmd = [
        "minimap2",
        "-x", "map-ont",
        "-t", str(threads),
        str(db_path),
        str(reads_path),
        "-o", str(paf_output_path)
    ]
    
    print(f"Aligning {reads_path} to {db_path}...")
    run_cmd(minimap_cmd)
    
    print(f"\nAlignment complete.")
    print(f"Output saved to: {paf_output_path}")
    print("Next step: Run scripts/4_generate_matrix.py to parse this file.")

if __name__ == "__main__":
    main()
