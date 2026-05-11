#!/usr/bin/env python3
"""
2b_trim_primers.py — Trim primer sequences from classified ONT reads.

Replaces the old Porechop ABI approach (which found ONT adapters, not primers)
with cutadapt using explicitly confirmed primer sequences.

Strategy per read:
  -g FWD_PRIMER  : find forward primer at 5'; drop ONT adapter chain + primer
  -a RC(REV)     : find RC of reverse primer at 3'; drop from there onward
  --revcomp      : try both orientations; output all reads as forward strand
  --error-rate   : 20% mismatches tolerated (accounts for ONT error rate ~10%)
  --overlap      : require ≥10 bp of primer match before trimming

Primer sequences are read from scripts/primers.tsv (tab-separated).
Edit that file to update or correct primers without touching this script.

Usage:
    python scripts/2b_trim_primers.py \\
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \\
        --markers 18S,COI \\
        --threads 8

    # Discard reads where no primer was found (stricter, fewer reads):
    python scripts/2b_trim_primers.py --input_dir ... --discard_untrimmed

    # Skip trimming entirely (pass-through):
    python scripts/2b_trim_primers.py --input_dir ... --skip
"""

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import find_tool

# Full IUPAC complement: handles degenerate bases Y, R, M, K, S, W, B, V, D, H, N
RC_MAP = str.maketrans(
    "ACGTRYMKSWBVDHNacgrtykmswbvdhn",
    "TGCAYRKMSWVBHDNtgcayrkmswvbhdn")


def rc(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(RC_MAP)[::-1]


# Built-in defaults — confirmed from k-mer analysis of actual reads.
# 18S Rev: identified from 3' end of forward-oriented reads (~63% hit rate).
# COI/JEDI: LCO1490 / HCO2198 (classic universal COI primers).
DEFAULT_PRIMERS = {
    "18S":  {"fwd": "ACCTGGTTGATCCTGCCAGT",       "rev": "TGTTACGACTTCACCTTCCTCTAAA"},
    "COI":  {"fwd": "GGTCAACAAATCATAAAGATATTGG",  "rev": "TAAACTTCAGGGTGACCAAAAAATCA"},
    "JEDI": {"fwd": "GTGYCAGCMGCCGCGGTAA",         "rev": "CCGYCAATTYMTTTRAGTTT"},    # 515F-Y / 926R
}


def load_primers(tsv_path: Path) -> dict:
    """
    Load primer sequences from primers.tsv.

    Falls back to DEFAULT_PRIMERS for any marker not listed in the file,
    so adding a new marker to the TSV is enough — no code change needed.
    """
    primers = dict(DEFAULT_PRIMERS)
    if not tsv_path.exists():
        return primers
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            marker = row.get("marker", "").strip()
            fwd    = row.get("fwd_primer", "").strip().upper()
            rev    = row.get("rev_primer", "").strip().upper()
            if marker and fwd and rev:
                primers[marker] = {"fwd": fwd, "rev": rev}
    return primers


def count_reads(path: Path) -> int:
    """Count reads in a gzipped FASTQ file (4 lines per read)."""
    n = 0
    try:
        with gzip.open(path, "rt") as fh:
            for _ in fh:
                n += 1
    except Exception:
        return 0
    return n // 4


def trim_barcode_marker(
    barcode_dir: Path,
    marker: str,
    fwd_primer: str,
    rev_primer: str,
    cutadapt_path: str,
    threads: int,
    discard_untrimmed: bool,
) -> tuple:
    """
    Trim one barcode/marker file with cutadapt.

    Backs up the original as .pretrim.fastq.gz before overwriting, so the
    QC notebook can compare pre- vs post-trim read lengths.
    The backup is created on first run, or refreshed when step 2 has produced
    newer classified reads (mtime guard prevents accidental overwrites).

    Returns (ok: bool, message: str).
    """
    fq_gz = barcode_dir / f"filtered_reads_{marker}.fastq.gz"
    if not fq_gz.exists():
        return False, "no input file"

    n_before = count_reads(fq_gz)
    pretrim  = barcode_dir / f"filtered_reads_{marker}.pretrim.fastq.gz"
    trimmed  = barcode_dir / f"filtered_reads_{marker}.trimmed.fastq.gz"

    # Refresh backup only when the classified reads are newer than the backup.
    if not pretrim.exists() or fq_gz.stat().st_mtime > pretrim.stat().st_mtime:
        shutil.copy2(fq_gz, pretrim)

    cmd = [
        cutadapt_path,
        "-g", fwd_primer,  # 5' adapter: everything before + including FWD primer
        "-a", rc(rev_primer),  # 3' adapter: RC(REV) as it appears at 3' of fwd reads
        "--revcomp",           # try RC of each read; output the better match as fwd
        "--error-rate", "0.2",
        "--overlap", "10",
        "--cores", str(threads),
        "--output", str(trimmed),
        str(fq_gz),
    ]
    if discard_untrimmed:
        cmd.insert(-1, "--discard-untrimmed")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0 or not trimmed.exists() or trimmed.stat().st_size == 0:
        trimmed.unlink(missing_ok=True)
        err = (result.stderr or result.stdout)[:300]
        return False, f"cutadapt failed (rc={result.returncode}): {err}"

    n_after = count_reads(trimmed)
    if n_after == 0:
        trimmed.unlink(missing_ok=True)
        return False, "all reads discarded — check primer sequences; original kept"

    shutil.move(str(trimmed), str(fq_gz))
    pct = 100 * n_after / max(n_before, 1)
    return True, f"{n_before} → {n_after} reads ({pct:.1f}% retained)"


def main():
    parser = argparse.ArgumentParser(
        description="Trim primers from classified reads using cutadapt.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing barcode subdirectories")
    parser.add_argument("--markers", default="18S,COI",
                        help="Comma-separated marker names (default: 18S,COI)")
    parser.add_argument("--threads", type=int, default=4,
                        help="CPU cores for cutadapt (default: 4)")
    parser.add_argument("--discard_untrimmed", action="store_true",
                        help="Discard reads where no primer was detected "
                             "(stricter; may lose coverage in sparse barcodes)")
    parser.add_argument("--skip", action="store_true",
                        help="Skip trimming entirely (pass-through)")
    args = parser.parse_args()

    if args.skip:
        print("Primer trimming skipped (--skip flag).")
        return 0

    cutadapt_path = find_tool("cutadapt") or find_tool("./env/bin/cutadapt")
    if not cutadapt_path:
        print("WARNING: cutadapt not found — skipping primer trimming.", file=sys.stderr)
        print("Install with:  env/bin/pip install cutadapt", file=sys.stderr)
        return 0

    print(f"cutadapt: {cutadapt_path}")

    primers_tsv = Path(__file__).resolve().parent / "primers.tsv"
    primers     = load_primers(primers_tsv)
    print(f"Primers: {primers_tsv if primers_tsv.exists() else 'built-in defaults'}")

    markers     = [m.strip() for m in args.markers.split(",")]
    input_dir   = Path(args.input_dir)

    barcode_dirs = sorted(
        d for d in input_dir.iterdir()
        if d.is_dir() and d.name.startswith("barcode")
    )
    if not barcode_dirs:
        print(f"No barcode directories found in {input_dir}")
        return 0

    for marker in markers:
        print(f"\n=== Marker: {marker} ===")

        if marker not in primers:
            print(f"  [SKIP] No primer sequences defined for '{marker}'.")
            print(f"         Add a row to {primers_tsv} or edit DEFAULT_PRIMERS.")
            continue

        fwd = primers[marker]["fwd"]
        rev = primers[marker]["rev"]
        print(f"  Fwd primer : {fwd}")
        print(f"  Rev primer : {rev}")
        print(f"  Rev RC (3'): {rc(rev)}")
        if args.discard_untrimmed:
            print("  --discard-untrimmed active: unmatched reads dropped")

        ok_count = skip_count = fail_count = 0
        for bdir in barcode_dirs:
            ok, msg = trim_barcode_marker(
                bdir, marker, fwd, rev, cutadapt_path,
                args.threads, args.discard_untrimmed,
            )
            tag = "[OK]  " if ok else ("[SKIP]" if "no input" in msg else "[WARN]")
            print(f"    {tag} {bdir.name}: {msg}")
            if ok:
                ok_count += 1
            elif "no input" in msg:
                skip_count += 1
            else:
                fail_count += 1

        print(f"  {marker}: {ok_count} trimmed, {skip_count} skipped, {fail_count} failed")

    return 0


if __name__ == "__main__":
    sys.exit(main())
