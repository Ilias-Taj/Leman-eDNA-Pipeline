#!/usr/bin/env python3
"""
2_classify_markers.py - CLASSIFY READS BY MARKER (18S vs COI)

Separates reads into 18S rRNA and COI based on sequence length only:
- 18S rRNA: 1,500-2,800 bp
- COI: 400-900 bp 
- Reads outside these ranges are discarded as noise

Prerequisite: Run 1_run_preprocessing.py first

For each barcode's filtered_reads.fastq.gz, outputs:
- filtered_reads_18S.fastq.gz
- filtered_reads_COI.fastq.gz
- marker_counts.txt

Usage:
    python scripts/2_classify_markers.py --input_dir out/Water_eDNA_18S_COI_14_01_26
"""

import argparse
import sys
import gzip
from pathlib import Path

# Sequence length thresholds (bp)
LENGTH_18S_MIN = 1500
LENGTH_18S_MAX = 2800
LENGTH_COI_MIN = 400   # Changed from 300 to avoid 18S contamination
LENGTH_COI_MAX = 900   # Changed from 1000 (more realistic for COI ~658bp)

def classify_read(seq, read_id):
    """
    Classify a read as 18S or COI based on sequence length only.
    Returns: "18S", "COI", or "ambiguous" (noise)
    """
    seq_len = len(seq)
    
    # 18S range: 1500-2800 bp
    if LENGTH_18S_MIN <= seq_len <= LENGTH_18S_MAX:
        return "18S"
    
    # COI range: 400-900 bp (stricter to avoid 18S contamination)
    elif LENGTH_COI_MIN <= seq_len <= LENGTH_COI_MAX:
        return "COI"
    
    # Reads outside both ranges are noise - discard them
    else:
        return "ambiguous"

def process_barcode(barcode_dir):
    """
    Process a single barcode directory.
    Splits filtered_reads.fastq.gz into 18S and COI files.
    """
    barcode_name = barcode_dir.name
    filtered_reads = barcode_dir / "filtered_reads.fastq.gz"
    
    if not filtered_reads.exists():
        print(f"  [SKIP] {barcode_name} - no filtered_reads.fastq.gz", file=sys.stderr)
        return None
    
    output_18s = barcode_dir / "filtered_reads_18S.fastq.gz"
    output_coi = barcode_dir / "filtered_reads_COI.fastq.gz"
    
    count_18s = 0
    count_coi = 0
    count_ambiguous = 0
    
    with gzip.open(filtered_reads, 'rt') as in_f:
        with gzip.open(output_18s, 'wt') as out_18s:
            with gzip.open(output_coi, 'wt') as out_coi:
                lines = []
                for line in in_f:
                    lines.append(line)
                    
                    # FASTQ: 4 lines per read
                    if len(lines) == 4:
                        header = lines[0]
                        seq = lines[1]
                        plus = lines[2]
                        qual = lines[3]
                        
                        marker = classify_read(seq.strip(), header.strip())
                        
                        if marker == "18S":
                            out_18s.write(header)
                            out_18s.write(seq)
                            out_18s.write(plus)
                            out_18s.write(qual)
                            count_18s += 1
                        elif marker == "COI":
                            out_coi.write(header)
                            out_coi.write(seq)
                            out_coi.write(plus)
                            out_coi.write(qual)
                            count_coi += 1
                        else:
                            count_ambiguous += 1
                        
                        lines = []
    
    total = count_18s + count_coi + count_ambiguous
    print(f"{barcode_name}: 18S={count_18s}, COI={count_coi}, ambiguous={count_ambiguous}, total={total}")
    
    # Write counts file
    counts_file = barcode_dir / "marker_counts.txt"
    with open(counts_file, 'w') as f:
        f.write(f"18S: {count_18s}\n")
        f.write(f"COI: {count_coi}\n")
        f.write(f"Ambiguous: {count_ambiguous}\n")
        f.write(f"Total: {total}\n")
    
    return {
        "barcode": barcode_name,
        "18S": count_18s,
        "COI": count_coi,
        "ambiguous": count_ambiguous,
        "total": total
    }

def main():
    parser = argparse.ArgumentParser(
        description="Classify eDNA reads by marker (18S vs COI)"
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Input directory containing barcode folders (e.g., out/Water_eDNA_18S_COI_14_01_26)"
    )
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"ERROR: Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    print("Classifying reads by marker (18S vs COI)...")
    print(f"Input directory: {input_dir}\n")
    
    results = []
    
    # Process each barcode directory
    for barcode_dir in sorted(input_dir.iterdir()):
        if not barcode_dir.is_dir():
            continue
        if barcode_dir.name in ("logs", "validation", "merged", "temp_clustering"):
            continue
        
        result = process_barcode(barcode_dir)
        if result:
            results.append(result)
    
    # Summary
    print("\n" + "="*60)
    print("MARKER CLASSIFICATION SUMMARY")
    print("="*60)
    
    total_18s = sum(r["18S"] for r in results)
    total_coi = sum(r["COI"] for r in results)
    total_ambiguous = sum(r["ambiguous"] for r in results)
    grand_total = sum(r["total"] for r in results)
    
    print(f"Total 18S reads: {total_18s} ({100*total_18s/grand_total:.1f}%)")
    print(f"Total COI reads: {total_coi} ({100*total_coi/grand_total:.1f}%)")
    print(f"Total ambiguous: {total_ambiguous} ({100*total_ambiguous/grand_total:.1f}%)")
    print(f"Grand total: {grand_total}")
    print("="*60)
    
    print("\nMarker-separated FASTQ files created for all barcodes.")
    print("Next: Run clustering separately for 18S and COI")

if __name__ == "__main__":
    main()
