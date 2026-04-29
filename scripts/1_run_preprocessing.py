#!/usr/bin/env python3
"""
1_run_preprocessing.py - QUALITY FILTERING

Filters reads by quality score using filtlong.
Handles multiple input FASTQ files (standard MinION output) by merging them.

Prerequisites:
    - Basecalled and demultiplexed FASTQ files from MinION
    - filtlong installed in environment

Usage:
    # Filter single or multiple files (shell glob)
    python scripts/1_run_preprocessing.py \
        --input_files data/Water_eDNA_18S_COI_14_01_26/fastq_pass/barcode01/*.fastq.gz \
        --output_dir out/Water_eDNA_18S_COI_14_01_26/barcode01 \
        --min_mean_q 20

Output:
    Creates filtered_reads.fastq.gz in the output directory

Note:
    Multiple input files are merged on-the-fly during filtering.
    No adapter/primer trimming is performed.
"""

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

# Ensure scripts/ is on the import path so utils.py can be found from any working directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import find_tool


def main():
    parser = argparse.ArgumentParser(
        description="Quality filter FASTQ reads using filtlong",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input_files", required=True, nargs='+',
                        help="Path(s) to input FASTQ or FASTQ.gz file(s)")
    parser.add_argument("--output_dir", required=True, 
                        help="Output directory (will be created if needed)")
    parser.add_argument("--min_length", type=int, default=0, 
                        help="Minimum read length (0=no length filter, default: 0)")
    parser.add_argument("--min_mean_q", type=int, default=20, 
                        help="Minimum mean quality score (default: 20)")
    parser.add_argument("--keep_percent", type=float, default=100.0, 
                        help="Keep top N%% of reads by quality (default: 100)")

    args = parser.parse_args()

    # Verify input files exist
    input_files = [Path(f) for f in args.input_files]
    for f in input_files:
        if not f.exists():
            print(f"ERROR: input file does not exist: {f}", file=sys.stderr)
            sys.exit(2)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    filtered_output = out_dir / "filtered_reads.fastq.gz"

    print("="*60)
    print("QUALITY FILTERING WITH FILTLONG")
    print("="*60)
    print(f"Input files:  {len(input_files)} file(s)")
    print(f"Output:       {filtered_output}")
    print(f"Min Q score:  {args.min_mean_q}")
    print(f"Min length:   {args.min_length if args.min_length > 0 else 'disabled'}")
    print(f"Keep percent: {args.keep_percent}%")
    print()

    # Locate filtlong executable (checks ./env/bin, system PATH, common locations)
    filtlong_path = find_tool("filtlong")

    # MinION outputs multiple FASTQ chunks per barcode.
    # Concatenate them into a single gzip stream before filtering.
    if len(input_files) > 1:
        print(f"Concatenating {len(input_files)} files...")
        temp_concat = out_dir / "temp_concatenated.fastq.gz"
        
        with gzip.open(temp_concat, 'wb') as outf:
            for f in input_files:
                with gzip.open(f, 'rb') as inf:
                    shutil.copyfileobj(inf, outf)
        
        input_for_filtlong = temp_concat
    else:
        input_for_filtlong = input_files[0]

    # Build filtlong command
    filt_cmd = [filtlong_path]
    if args.min_length > 0:
        filt_cmd += ["--min_length", str(args.min_length)]
    if args.min_mean_q > 0:
        filt_cmd += ["--min_mean_q", str(args.min_mean_q)]
    if args.keep_percent < 100:
        filt_cmd += ["--keep_percent", str(args.keep_percent)]
    filt_cmd.append(str(input_for_filtlong))

    print("Running filtlong...")

    try:
        with gzip.open(filtered_output, "wb") as gzout:
            # Run filtlong
            proc = subprocess.Popen(filt_cmd, stdout=subprocess.PIPE)
            shutil.copyfileobj(proc.stdout, gzout)
            proc.stdout.close()
            ret = proc.wait()
            
            if ret != 0:
                print(f"ERROR: filtlong failed with code {ret}", file=sys.stderr)
                sys.exit(ret)

    except Exception as e:
        print(f"ERROR during execution: {e}", file=sys.stderr)
        sys.exit(1)
    
    finally:
        # Clean up temp file if we created one
        if len(input_files) > 1 and temp_concat.exists():
            temp_concat.unlink()

    print(f"\n[OK] Filtered reads saved to: {filtered_output}")
    print("\nNext step: Classify markers with scripts/2_classify_markers.py")


if __name__ == "__main__":
    main()
