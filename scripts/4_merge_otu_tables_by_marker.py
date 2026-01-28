#!/usr/bin/env python3
"""
4_merge_otu_tables_by_marker.py - GENERATE SEPARATE ABUNDANCE MATRICES PER MARKER

Creates abundance matrices separately for 18S and COI markers.
Each matrix is in OTU × sample format with both raw and relative abundance.

Prerequisite: Run 3_run_clustering_by_marker.py first

Usage:
    python scripts/4_merge_otu_tables_by_marker.py --input_dir out/Water_eDNA_18S_COI_14_01_26

Output:
    Saves 4 CSV files to {input_dir}/merged/:
    - otu_abundance_matrix_18S.csv (raw counts)
    - otu_relative_abundance_18S.csv (proportions 0-1)
    - otu_abundance_matrix_COI.csv (raw counts)
    - otu_relative_abundance_COI.csv (proportions 0-1)
"""

import sys
import argparse
from pathlib import Path
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description="Generate abundance matrices from OTU assignments")
    parser.add_argument("--input_dir", required=True, help="Path to the run directory (e.g. out/Run_Name)")
    args = parser.parse_args()
    
    run_dir = Path(args.input_dir)
    
    if not run_dir.exists():
        print(f"ERROR: Directory not found: {run_dir}", file=sys.stderr)
        sys.exit(1)
    
    merged_dir = run_dir / "merged"
    merged_dir.mkdir(exist_ok=True)
    
    # Load global assignment file
    assignment_file = run_dir / "global_otu_assignment.txt"
    if not assignment_file.exists():
        print(f"ERROR: Assignment file not found: {assignment_file}", file=sys.stderr)
        print("Did you run Script 3 (Clustering)?", file=sys.stderr)
        sys.exit(1)
    
    print(f"Reading assignments from: {assignment_file}")
    
    # Structure: {marker: {otu: {barcode: count}}}
    otu_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    all_barcodes = set()
    
    # Robust reading: handle potential header or empty lines
    try:
        with open(assignment_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("read_name"):  # Skip empty lines and header
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                
                # columns: read_name, otu_id, barcode, marker
                otu_id, barcode, marker = parts[1], parts[2], parts[3]
                
                otu_counts[marker][otu_id][barcode] += 1
                all_barcodes.add(barcode)
    except Exception as e:
        print(f"ERROR reading assignment file: {e}", file=sys.stderr)
        sys.exit(1)
    
    all_barcodes = sorted(list(all_barcodes))
    print(f"Found {len(all_barcodes)} samples: {all_barcodes}")
    
    # Generate separate matrices for each marker
    # We iterate specifically over 18S and COI to enforce separation
    for marker in ["18S", "COI"]:
        if marker not in otu_counts:
            print(f"\n[INFO] No data found for marker: {marker} (skipping)")
            continue
            
        print(f"\nGenerating {marker} matrices...")
        
        # Get all OTUs for this marker
        otus = sorted(otu_counts[marker].keys())
        print(f"  > {len(otus)} OTUs found")
        
        # 1. Write Raw Abundance Matrix
        matrix_file = merged_dir / f"otu_abundance_matrix_{marker}.csv"
        with open(matrix_file, 'w') as f:
            # Header: OTU_ID,barcode01,barcode02...
            f.write("OTU_ID," + ",".join(all_barcodes) + "\n")
            
            for otu in otus:
                # List comprehension ensures columns align with header
                counts = [str(otu_counts[marker][otu].get(bc, 0)) for bc in all_barcodes]
                f.write(f"{otu}," + ",".join(counts) + "\n")
        print(f"  ✓ Saved raw counts: {matrix_file.name}")
        
        # 2. Write Relative Abundance Matrix
        rel_matrix_file = merged_dir / f"otu_relative_abundance_{marker}.csv"
        
        # Calculate column totals (depth per sample per marker)
        # Crucial: Normalization must be done per sample!
        sample_totals = {}
        for bc in all_barcodes:
            total = sum(otu_counts[marker][otu].get(bc, 0) for otu in otus)
            sample_totals[bc] = total

        with open(rel_matrix_file, 'w') as f:
            f.write("OTU_ID," + ",".join(all_barcodes) + "\n")
            
            for otu in otus:
                rel_counts = []
                for bc in all_barcodes:
                    count = otu_counts[marker][otu].get(bc, 0)
                    total = sample_totals[bc]
                    
                    # Robust division: avoid ZeroDivisionError if a sample has 0 reads for this marker
                    if total > 0:
                        norm_val = count / total
                    else:
                        norm_val = 0.0
                    rel_counts.append(f"{norm_val:.6f}")
                
                f.write(f"{otu}," + ",".join(rel_counts) + "\n")
        print(f"  ✓ Saved relative abundance: {rel_matrix_file.name}")
    
    print(f"\nSuccess. Matrices saved in: {merged_dir}")

if __name__ == "__main__":
    main()
