#!/usr/bin/env python3
"""
2_run_preprocessing.py

A robust preprocessing pipeline for MinION eDNA metabarcoding data.
This script is based on common literature workflows.

Steps:
 - Trim adapters/primers with cutadapt (single-end amplicon trimming)
 - Filter trimmed reads with filtlong

Usage example:
    python 2_run_preprocessing.py \
        --input_file test_reads.fastq \
        --output_dir out \
        --primer_f GTCGGTAAAACTCGTGCCAGC \
        --primer_r CATAGENTGACGTGGGRGGT \
        --min_length 100 \
        --keep_percent 90
"""

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

def run_cmd(cmd, shell=False):
    """Prints and runs a command, exiting on failure."""
    if shell:
        print(f"Running shell command: {cmd}")
    else:
        print("Running:", " ".join(map(str, cmd)))

    try:
        subprocess.run(cmd, check=True, shell=shell, text=True, executable='/bin/bash')
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)
    except FileNotFoundError as e:
        tool = cmd.split()[0]
        print(f"ERROR: Command not found: {tool}", file=sys.stderr)
        print(f"Please make sure '{tool}' is installed and in your PATH.", file=sys.stderr)
        sys.exit(1)

def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S',
                  'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                  'D': 'H', 'H': 'D', 'V': 'B'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def main():
    parser = argparse.ArgumentParser(description="Preprocessing: cutadapt -> Filtlong")
    parser.add_argument("--input_file", required=True, help="Path to input FASTQ file (can be .gz)")
    parser.add_argument("--output_dir", required=True, help="Directory to save processed files")
    parser.add_argument("--primer_f", required=True, help="Forward primer sequence (5' to 3')")
    parser.add_argument("--primer_r", required=True, help="Reverse primer sequence (5' to 3')")
    parser.add_argument("--min_length", type=int, default=100, help="Minimum read length to keep (post-trimming)")
    parser.add_argument("--keep_percent", type=float, default=90.0, help="Percent of best reads to keep (e.g., 90.0)")
    parser.add_argument("--threads", type=int, default=4, help="Number of CPU threads to use")

    args = parser.parse_args()

    input_path = Path(args.input_file)
    out_dir = Path(args.output_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    final_output_path = out_dir / "filtered_reads.fastq.gz"
    cutadapt_report = out_dir / "cutadapt_report.txt"
    trimmed_fastq = out_dir / "trimmed.fastq"
    trimmed_fasta = out_dir / "trimmed.fasta"

    fwd_primer = args.primer_f
    rev_primer = args.primer_r
    fwd_primer_rc = get_reverse_complement(fwd_primer)
    rev_primer_rc = get_reverse_complement(rev_primer)

    print("--- Primer Setup ---")
    print(f"Fwd Primer (5'-3'):       {fwd_primer}")
    print(f"Rev Primer (5'-3'):       {rev_primer}")
    print(f"Fwd Primer (RC):         {fwd_primer_rc}")
    print(f"Rev Primer (RC):         {rev_primer_rc}")
    print("----------------------")

    print("\n--- Starting Preprocessing Pipeline ---")

    def run_cutadapt(output_path):
        cmd = [
            "cutadapt",
            "-j", str(args.threads),
            "--discard-untrimmed",
            "--minimum-length", str(args.min_length),
            "-g", fwd_primer,
            "-a", rev_primer_rc,
            "-g", rev_primer,
            "-a", fwd_primer_rc,
            "-o", str(output_path),
            str(input_path)
        ]
        print("Running cutadapt:", " ".join(cmd))
        with open(cutadapt_report, "w") as rpt:
            try:
                subprocess.run(cmd, check=True, stderr=rpt)
            except subprocess.CalledProcessError as e:
                print(f"ERROR: cutadapt failed with code {e.returncode}", file=sys.stderr)
                sys.exit(e.returncode)

    is_fasta = False
    with open(str(input_path), 'rb') as fh:
        start = fh.read(1)
        if start == b'>':
            is_fasta = True

    if is_fasta:
        run_cutadapt(trimmed_fasta)
        import gzip
        print(f"Converting trimmed FASTA to FASTQ (dummy qualities) -> {trimmed_fastq}")
        with open(trimmed_fasta, "rt") as inf, open(trimmed_fastq, "wt") as outf:
            name = None
            seq_lines = []
            for line in inf:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if name is not None:
                        seq = ''.join(seq_lines)
                        outf.write(f"@{name}\n")
                        outf.write(seq + "\n")
                        outf.write("+\n")
                        outf.write('I' * len(seq) + "\n")
                    name = line[1:].split()[0]
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if name is not None:
                seq = ''.join(seq_lines)
                outf.write(f"@{name}\n")
                outf.write(seq + "\n")
                outf.write("+\n")
                outf.write('I' * len(seq) + "\n")

        import shutil
        with open(trimmed_fastq, 'rb') as f_in, gzip.open(str(final_output_path), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        print(f"Wrote gzipped FASTQ to: {final_output_path}")

    else:
        run_cutadapt(trimmed_fastq)
        filtlong_cmd = [
            "filtlong",
            "--min_length", str(args.min_length),
            "--keep_percent", str(args.keep_percent),
            str(trimmed_fastq)
        ]

        print("Running filtlong:", " ".join(filtlong_cmd))
        import gzip
        with gzip.open(str(final_output_path), 'wb') as gzout:
            p = subprocess.Popen(filtlong_cmd, stdout=subprocess.PIPE)
            try:
                import shutil
                shutil.copyfileobj(p.stdout, gzout)
            finally:
                p.stdout.close()
                ret = p.wait()
                if ret != 0:
                    print(f"ERROR: filtlong exited with code {ret}", file=sys.stderr)
                    sys.exit(ret)

        print(f"Wrote filtered (gzipped) FASTQ to: {final_output_path}")

    print("\n--- Preprocessing Complete ---")
    print(f"Cutadapt report saved to: {cutadapt_report}")
    print(f"Final processed reads saved to: {final_output_path}")
    print("This file is now ready for alignment with '3_run_alignment.py'.")

if __name__ == "__main__":
    main()
