#!/usr/bin/env python3
"""
merge_and_compare.py

Merge per-sample species_counts and validation reports, compare top MinION hit to Sanger reference
found in each sample's SANGER/*.fasta file.

Outputs:
 - out/merged/species_by_sample.csv  (reference x sample relative abundance matrix)
 - out/merged/comparison_report.csv  (one row per sample: sanger_id, top_minion_ref, match, validation_status, identity)

"""
import csv
from pathlib import Path
import sys

ROOT = Path("raw_reads")
OUT = Path("out")
OUT.mkdir(parents=True, exist_ok=True)

samples = [p.name for p in sorted(ROOT.iterdir()) if p.is_dir()]

# Collect all reference ids found across samples
refs = set()
per_sample_counts = {}

for s in samples:
    paf_counts = OUT / s / "species_counts.csv"
    counts = {}
    if paf_counts.exists():
        with open(paf_counts, "rt") as fh:
            rdr = csv.reader(fh)
            hdr = next(rdr)
            for row in rdr:
                if not row:
                    continue
                ref = row[0]
                try:
                    c = int(row[1])
                except Exception:
                    c = 0
                counts[ref] = c
                refs.add(ref)
    per_sample_counts[s] = counts

# Build matrix (rows=ref, cols=samples) with relative abundance
merged_dir = OUT / "merged"
merged_dir.mkdir(parents=True, exist_ok=True)
matrix_file = merged_dir / "species_by_sample.csv"

all_refs = sorted(refs)

with open(matrix_file, "wt", newline="") as mf:
    w = csv.writer(mf)
    header = ["reference_id"] + samples
    w.writerow(header)
    for r in all_refs:
        row = [r]
        for s in samples:
            total = sum(per_sample_counts.get(s, {}).values())
            c = per_sample_counts.get(s, {}).get(r, 0)
            rel = c / total if total > 0 else 0.0
            row.append(f"{rel:.6f}")
        w.writerow(row)

# Now build comparison report
comp_file = merged_dir / "comparison_report.csv"
with open(comp_file, "wt", newline="") as cf:
    w = csv.writer(cf)
    w.writerow(["sample","sanger_id","top_minion_ref","top_count","top_rel_abundance","minion_top_matches_sanger","validation_status","validation_identity_pct"])

    for s in samples:
        # sanger id: read first header token from SANGER fasta in sample folder
        sanger_dir = ROOT / s / "SANGER"
        sanger_id = ""
        if sanger_dir.exists():
            fa_files = list(sanger_dir.glob("*.fasta"))
            if fa_files:
                # read first fasta header
                with open(fa_files[0], "rt") as fh:
                    for line in fh:
                        if line.startswith(">"):
                            sanger_id = line[1:].strip().split()[0]
                            break

        # top minion ref from species_counts
        counts = per_sample_counts.get(s, {})
        top_ref = ""
        top_count = 0
        top_rel = 0.0
        total = sum(counts.values())
        if counts:
            top_ref, top_count = max(counts.items(), key=lambda x: x[1])
            top_rel = top_count / total if total > 0 else 0.0

        # read validation status file if exists
        val_file = OUT / "validation" / s / "validation_report.csv"
        val_status = ""
        val_identity = ""
        if val_file.exists():
            with open(val_file, "rt") as vf:
                rdr = csv.reader(vf)
                hdr = next(rdr, None)
                for row in rdr:
                    if not row:
                        continue
                    # choose row corresponding to top_ref if present, else first row
                    if row[0] == top_ref or top_ref == "":
                        val_status = row[5] if len(row) > 5 else ""
                        val_identity = row[4] if len(row) > 4 else ""
                        break

        # matching logic: consider a match if exact token match or prefix match
        match = False
        if sanger_id and top_ref:
            if top_ref == sanger_id:
                match = True
            elif top_ref.startswith(sanger_id):
                match = True
            elif sanger_id.startswith(top_ref):
                match = True
            else:
                # sometimes species name is part of ref id before a '|'
                if '|' in top_ref:
                    if top_ref.split('|',1)[0] == sanger_id:
                        match = True

        w.writerow([s, sanger_id, top_ref, top_count, f"{top_rel:.6f}", str(match), val_status, val_identity])

print("Merged matrix:", matrix_file)
print("Comparison report:", comp_file)
