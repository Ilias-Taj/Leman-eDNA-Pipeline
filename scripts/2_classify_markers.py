#!/usr/bin/env python3
"""
2_classify_markers.py - CLASSIFY READS BY MARKER

Separates reads into marker bins based on sequence length:
- 18S rRNA: 1,500-2,800 bp
- COI:      500-900 bp   (standard Folmer/Leray primers, ~658 bp amplicon)
- JEDI:     250-500 bp   (JEDI COI primers, ~460 bp amplicon)

Use --markers to select which markers to search for (default: 18S,COI).
Reads outside all selected marker ranges are discarded as noise.

Prerequisite: Run 1_run_preprocessing.py first

For each barcode's filtered_reads.fastq.gz, outputs:
- filtered_reads_{marker}.fastq.gz  (one per selected marker)
- marker_counts.txt

Usage:
    # Water data (18S + COI)
    python scripts/2_classify_markers.py --input_dir out/Water_eDNA_18S_COI_14_01_26 --markers 18S,COI

    # Soil data (JEDI + COI)
    python scripts/2_classify_markers.py --input_dir out/Soil_eDNA_JEDI_COI_14_01_26 --markers JEDI,COI

    # All three markers
    python scripts/2_classify_markers.py --input_dir out/run_name --markers 18S,COI,JEDI
"""

import argparse
import sys
import gzip
from pathlib import Path

# Sequence length thresholds (bp) for each marker
# Derived from empirical sequencing data:
#   Water run N50=1.6kb (18S dominates), Soil run N50=636bp (JEDI+COI mix)
MARKER_RANGES = {
    "18S":  (1500, 2800),   # 18S rRNA (~1.8 kb amplicon)
    "COI":  (500,  900),    # Standard COI Folmer/Leray (~658 bp amplicon)
    "JEDI": (250,  500),    # JEDI short COI primers (~460 bp amplicon)
}

VALID_MARKERS = list(MARKER_RANGES.keys())

def classify_read(seq, read_id, active_markers):
    """
    Classify a read based on sequence length into one of the active markers.
    Returns the marker name or "ambiguous" if it doesn't fit any range.
    """
    seq_len = len(seq)
    
    for marker in active_markers:
        lo, hi = MARKER_RANGES[marker]
        if lo <= seq_len <= hi:
            return marker
    
    return "ambiguous"

def process_barcode(barcode_dir, active_markers):
    """
    Process a single barcode directory.
    Splits filtered_reads.fastq.gz into per-marker files.
    """
    barcode_name = barcode_dir.name
    filtered_reads = barcode_dir / "filtered_reads.fastq.gz"
    
    if not filtered_reads.exists():
        print(f"  [SKIP] {barcode_name} - no filtered_reads.fastq.gz", file=sys.stderr)
        return None
    
    # Open output files for each active marker
    output_handles = {}
    for marker in active_markers:
        output_handles[marker] = gzip.open(barcode_dir / f"filtered_reads_{marker}.fastq.gz", 'wt')
    
    counts = {m: 0 for m in active_markers}
    counts["ambiguous"] = 0
    
    with gzip.open(filtered_reads, 'rt') as in_f:
        lines = []
        for line in in_f:
            lines.append(line)
            
            # FASTQ: 4 lines per read
            if len(lines) == 4:
                header = lines[0]
                seq = lines[1]
                plus = lines[2]
                qual = lines[3]
                
                marker = classify_read(seq.strip(), header.strip(), active_markers)
                
                if marker in output_handles:
                    output_handles[marker].write(header)
                    output_handles[marker].write(seq)
                    output_handles[marker].write(plus)
                    output_handles[marker].write(qual)
                    counts[marker] += 1
                else:
                    counts["ambiguous"] += 1
                
                lines = []
    
    # Close all output files
    for fh in output_handles.values():
        fh.close()
    
    total = sum(counts.values())
    parts = ", ".join(f"{m}={counts[m]}" for m in active_markers)
    print(f"{barcode_name}: {parts}, ambiguous={counts['ambiguous']}, total={total}")
    
    # Write counts file
    counts_file = barcode_dir / "marker_counts.txt"
    with open(counts_file, 'w') as f:
        for m in active_markers:
            f.write(f"{m}: {counts[m]}\n")
        f.write(f"Ambiguous: {counts['ambiguous']}\n")
        f.write(f"Total: {total}\n")
    
    return {
        "barcode": barcode_name,
        **{m: counts[m] for m in active_markers},
        "ambiguous": counts["ambiguous"],
        "total": total
    }

def main():
    parser = argparse.ArgumentParser(
        description="Classify eDNA reads by marker based on amplicon length"
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Input directory containing barcode folders (e.g., out/Water_eDNA_18S_COI_14_01_26)"
    )
    parser.add_argument(
        "--markers",
        default="18S,COI",
        help=(
            "Comma-separated list of markers to classify. "
            f"Valid markers: {', '.join(VALID_MARKERS)}. "
            "Default: 18S,COI. For soil JEDI data use: JEDI,COI"
        )
    )
    
    args = parser.parse_args()
    
    # Parse and validate markers
    active_markers = [m.strip().upper() for m in args.markers.split(",")]
    for m in active_markers:
        if m not in VALID_MARKERS:
            print(f"ERROR: Unknown marker '{m}'. Valid markers: {', '.join(VALID_MARKERS)}", file=sys.stderr)
            sys.exit(1)
    
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"ERROR: Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Classifying reads by marker ({', '.join(active_markers)})...")
    print(f"Input directory: {input_dir}")
    print(f"Length ranges:")
    for m in active_markers:
        lo, hi = MARKER_RANGES[m]
        print(f"  {m}: {lo}-{hi} bp")
    print()
    
    results = []
    
    # Process each barcode directory
    for barcode_dir in sorted(input_dir.iterdir()):
        if not barcode_dir.is_dir():
            continue
        if barcode_dir.name in ("logs", "validation", "merged", "temp_clustering", "taxonomy", "taxonomy_summary", "blast_results"):
            continue
        
        result = process_barcode(barcode_dir, active_markers)
        if result:
            results.append(result)
    
    # Summary
    print("\n" + "="*60)
    print("MARKER CLASSIFICATION SUMMARY")
    print("="*60)
    
    grand_total = sum(r["total"] for r in results)
    if grand_total == 0:
        print("No reads processed!")
        return
    
    for m in active_markers:
        total_m = sum(r.get(m, 0) for r in results)
        print(f"Total {m} reads: {total_m} ({100*total_m/grand_total:.1f}%)")
    total_ambiguous = sum(r["ambiguous"] for r in results)
    print(f"Total ambiguous: {total_ambiguous} ({100*total_ambiguous/grand_total:.1f}%)")
    print(f"Grand total: {grand_total}")
    print("="*60)
    
    print(f"\nMarker-separated FASTQ files created for all barcodes.")
    print(f"Active markers: {', '.join(active_markers)}")
    print("Next: Run clustering separately for each marker")

if __name__ == "__main__":
    main()
