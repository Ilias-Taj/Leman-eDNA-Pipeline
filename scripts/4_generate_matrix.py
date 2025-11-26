#!/usr/bin/env python3
"""
4_generate_matrix.py

Parse a PAF produced by minimap2 and produce a simple species-by-sample CSV.

Usage:
  python scripts/4_generate_matrix.py --paf out/alignments.paf --output out/species_counts.csv

Output columns:
  reference_id,read_count,relative_abundance,presence

This is a minimal parser: it counts best hits per read by choosing the highest "matching bases" value (PAF col 10).
"""

import argparse
from collections import defaultdict
from pathlib import Path
import csv


def parse_paf(paf_path):
    best = {}
    with open(paf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            qname = cols[0]
            try:
                n_match = int(cols[9])
            except ValueError:
                n_match = 0
            tname = cols[5]

            if qname not in best or n_match > best[qname][0]:
                best[qname] = (n_match, tname)

    counts = defaultdict(int)
    total = 0
    for q, (n, t) in best.items():
        counts[t] += 1
        total += 1

    return counts, total


def main():
    parser = argparse.ArgumentParser(description="Parse PAF and create species counts CSV")
    parser.add_argument("--paf", required=True, help="PAF file produced by minimap2")
    parser.add_argument("--output", required=True, help="CSV output path")
    args = parser.parse_args()

    paf = Path(args.paf)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    counts, total = parse_paf(paf)

    with open(out, "w", newline="") as csvf:
        w = csv.writer(csvf)
        w.writerow(["reference_id", "read_count", "relative_abundance", "presence"])
        for ref, c in sorted(counts.items(), key=lambda x: -x[1]):
            rel = c / total if total > 0 else 0.0
            pres = 1 if c > 0 else 0
            w.writerow([ref, c, f"{rel:.6f}", pres])

    print(f"Wrote species counts for {total} reads to: {out}")


if __name__ == "__main__":
    main()
