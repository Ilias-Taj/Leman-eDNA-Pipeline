#!/usr/bin/env python3
"""
Draft basecalling step using Dorado and .pod5 raw signals.

This script is a draft. It locates .pod5 files for a sample and calls Dorado
to produce FASTQ output. The command builder is a placeholder and should be
adjusted to match the installed Dorado CLI and model names/paths on the
target system.

Example:
    python scripts/0_dorado_basecall.py --pod5_dir raw_reads/<sample>/pod5 --output out/<sample>/raw.fastq.gz --model dna_r10.4
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd, dry_run=False):
    # Print the command and execute it unless in dry-run mode.
    print("DRY RUN:" if dry_run else "Running:", " ".join(map(str, cmd)))
    if dry_run:
        # Do not execute in dry-run mode (is that what we want?)
        return 0
    try:
        # Execute the command and raise on non-zero exit.
        return subprocess.run(cmd, check=True).returncode
    except subprocess.CalledProcessError as e:
        print(f"ERROR: dorado/basecalling command failed: {e}", file=sys.stderr)
        return e.returncode


def find_pod5_files(pod5_dir):
    # Return a sorted list of .pod5 files (searches recursively).
    p = Path(pod5_dir)
    if not p.exists():
        raise FileNotFoundError(f"pod5 directory not found: {pod5_dir}")
    files = sorted(p.rglob("*.pod5"))
    return files


def build_dorado_command(pod5_input, out_fastq, model, device=None, threads=4, extra_args=None):
    # Build a placeholder dorado command. Adjust to match local dorado CLI.
    cmd = ["dorado", "basecaller"]
    if model:
        cmd += ["-m", str(model)]
    if device:
        cmd += ["--device", str(device)]
    # use --workers as a common name.
    cmd += ["--workers", str(threads)]
    cmd += [str(pod5_input)]
    cmd += ["-o", str(out_fastq)]
    if extra_args:
        cmd += extra_args
    return cmd


def main():
    parser = argparse.ArgumentParser(description="Draft: basecall .pod5 files with Dorado (commented)")
    parser.add_argument("--pod5_dir", required=True, help="Directory containing .pod5 files for a sample (or a single .pod5 file)")
    parser.add_argument("--output", required=True, help="Output FASTQ path (can be .gz)")
    parser.add_argument("--model", required=False, default=None, help="Dorado model name or path (e.g., dna_r10.4)")
    parser.add_argument("--device", required=False, default=None, help="Device string (e.g., cuda:0 or cpu)")
    parser.add_argument("--threads", type=int, default=4, help="Number of workers/threads for dorado")
    parser.add_argument("--dry_run", action="store_true", help="Print commands but do not execute them")
    args = parser.parse_args()

    pod5 = Path(args.pod5_dir)
    outp = Path(args.output)

    # Create parent dirs for the requested output.
    outp.parent.mkdir(parents=True, exist_ok=True)

    # Find .pod5 files in the provided directory.
    files = find_pod5_files(pod5)
    if not files:
        print(f"No .pod5 files found in: {pod5}", file=sys.stderr)
        sys.exit(2)

    # Use the pod5 directory path as dorado input (may need adjustment per install).
    pod5_input = str(pod5)

    # Build the dorado command / adjust to local install as needed.
    cmd = build_dorado_command(pod5_input, outp, args.model, device=args.device, threads=args.threads)

    # Print debugging info and the command to be run.
    print("Found pod5 files (example):", files[:3])
    print("Command to run (may need manual edit):")
    print(" ".join(map(str, cmd)))

    # Run (or dry-run)
    rc = run_cmd(cmd, dry_run=args.dry_run)
    if rc != 0:
        print("Basecalling command returned non-zero exit code.", file=sys.stderr)
        sys.exit(rc)

    # Check that the output file exists and report size.
    if outp.exists():
        sz = outp.stat().st_size
        print(f"Wrote basecalled FASTQ to: {outp} (size: {sz} bytes)")
    else:
        print(f"Expected output not found: {outp}", file=sys.stderr)
        sys.exit(3)


if __name__ == "__main__":
    main()
