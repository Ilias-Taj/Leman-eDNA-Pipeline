#!/usr/bin/env python3
"""
2b_trim_primers.py — Primer-based reclassification + trimming of ONT reads.

Runs after step 2 (isONclust length classification) and replaces it as the
source of truth for marker assignment.

PHASE 1 — Reclassify (primer scan):
  Pool all reads from all marker files for each barcode.
  Scan the first 200bp (and RC of first 200bp) of each read for every
  primer sequence in primers.tsv, with up to 20% mismatch tolerance.
  Assign the read to the marker whose primer is found.
  Discard reads where no primer is detected (not target amplicons:
  likely chimeras, cross-contamination, or non-target DNA).

PHASE 2 — Trim (cutadapt):
  -g FWD_PRIMER  : remove everything up to and including the forward primer
  -a RC(REV)     : remove RC of reverse primer and everything after it
  --revcomp      : normalise all reads to forward orientation

Primer sequences are read from scripts/primers.tsv (tab-separated):
  marker   fwd_primer   rev_primer

Edit primers.tsv to update primers without changing this script.

Usage:
    python scripts/2b_trim_primers.py \
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \
        --markers 18S,COI \
        --threads 8

    # Keep reads even when no primer found (no reclassification, just trim):
    python scripts/2b_trim_primers.py --input_dir ... --no_reclassify

    # Skip trimming entirely (pass-through):
    python scripts/2b_trim_primers.py --input_dir ... --skip
"""

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import find_tool

# Full IUPAC complement table for rc()
_IUPAC_COMP = str.maketrans(
    "ACGTRYMKSWBVDHNacgtryMkswbvdhn",
    "TGCAYRKMSWVBHDNtgcayrMkswvbhdn")


def rc(seq: str) -> str:
    """Reverse complement of a DNA sequence (handles IUPAC degenerate bases)."""
    return seq.translate(_IUPAC_COMP)[::-1]


def iupac_to_regex(primer: str) -> re.Pattern:
    """
    Expand an IUPAC degenerate primer to a regex pattern.
    Returns a compiled pattern for the primer alone (no leading/trailing context).
    """
    _map = {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "[AG]", "Y": "[CT]", "M": "[AC]", "K": "[GT]",
        "S": "[CG]", "W": "[AT]", "B": "[CGT]", "V": "[ACG]",
        "D": "[AGT]", "H": "[ACT]", "N": "[ACGT]",
    }
    return re.compile("".join(_map.get(b.upper(), b) for b in primer))


def primer_in_region(region: str, pattern: re.Pattern, max_err_rate: float = 0.20) -> bool:
    """
    Check whether a primer pattern appears in *region* with up to max_err_rate
    mismatches.  Uses exact regex first (fast), then a sliding-window Hamming
    distance check for mismatches if the exact match fails.
    """
    if pattern.search(region):
        return True

    primer_str = pattern.pattern  # the raw (expanded) string — we use original
    # Fall back to Hamming distance on each window of the region
    # We need the original primer string (not the regex expansion) for this.
    # We'll get it from the caller via a closure — see classify_read().
    return False


def hamming_match(region: str, primer: str, max_err_rate: float = 0.20) -> bool:
    """
    Sliding-window Hamming distance match (IUPAC-aware on primer side).
    Returns True if primer occurs anywhere in region with ≤ max_err_rate
    mismatches.  IUPAC codes in the primer always match any of their bases.
    """
    plen = len(primer)
    max_err = int(plen * max_err_rate)
    primer_up = primer.upper()

    _iupac_sets = {
        "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
        "R": {"A","G"}, "Y": {"C","T"}, "M": {"A","C"}, "K": {"G","T"},
        "S": {"C","G"}, "W": {"A","T"}, "B": {"C","G","T"}, "V": {"A","C","G"},
        "D": {"A","G","T"}, "H": {"A","C","T"}, "N": {"A","C","G","T"},
    }

    region_up = region.upper()
    for start in range(len(region_up) - plen + 1):
        window = region_up[start:start + plen]
        mismatches = 0
        for rb, pb in zip(window, primer_up):
            allowed = _iupac_sets.get(pb, {pb})
            if rb not in allowed:
                mismatches += 1
                if mismatches > max_err:
                    break
        if mismatches <= max_err:
            return True
    return False


def classify_read(seq: str, primers: dict, search_len: int = 200,
                  max_err_rate: float = 0.20) -> str | None:
    """
    Determine which marker (if any) owns this read by searching for primer
    sequences in the first *search_len* bp and the RC of the last *search_len* bp.

    Returns the marker name (e.g. '18S') or None if no primer found.
    """
    fwd_region = seq[:search_len].upper()
    rev_region = rc(seq[-search_len:]).upper()  # as if read were in fwd orientation

    for marker, info in primers.items():
        fwd_primer = info["fwd"]
        # Check forward primer in either region (some reads start with fwd, some end)
        if hamming_match(fwd_region, fwd_primer, max_err_rate):
            return marker
        if hamming_match(rev_region, fwd_primer, max_err_rate):
            return marker
    return None


DEFAULT_PRIMERS = {
    "18S":  {"fwd": "ACCTGGTTGATCCTGCCAGT",       "rev": "TGTTACGACTTCACCTTCCTCTAAA"},
    "COI":  {"fwd": "GGTCAACAAATCATAAAGATATTGG",  "rev": "TAAACTTCAGGGTGACCAAAAATCA"},
    "JEDI": {"fwd": "GTGYCAGCMGCCGCGGTAA",         "rev": "CCGYCAATTYMTTTRAGTTT"},
}


def load_primers(tsv_path: Path) -> dict:
    """Load primer sequences from primers.tsv, falling back to defaults."""
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


def read_fastq_gz(path: Path) -> list:
    """Return list of (header, seq, plus, qual) tuples from a gzipped FASTQ."""
    records = []
    try:
        with gzip.open(path, "rt") as f:
            lines = f.readlines()
        for i in range(0, len(lines) - 3, 4):
            records.append((lines[i].rstrip("\n"),
                            lines[i+1].rstrip("\n"),
                            lines[i+2].rstrip("\n"),
                            lines[i+3].rstrip("\n")))
    except Exception:
        pass
    return records


def write_fastq_gz(records: list, path: Path) -> None:
    """Write list of (header, seq, plus, qual) tuples to a gzipped FASTQ."""
    with gzip.open(path, "wt") as f:
        for h, s, p, q in records:
            f.write(f"{h}\n{s}\n{p}\n{q}\n")


def reclassify_barcode(barcode_dir: Path, markers: list, primers: dict,
                       max_err_rate: float = 0.20) -> dict:
    """
    Pool all reads from all marker files in this barcode directory.
    Scan each read for primer sequences and assign to the correct marker.
    Overwrite all marker files with the reclassified reads.
    Discards reads where no primer is detected.

    Returns a stats dict: {marker: {"kept": int, "discarded": int}, "total": int}
    """
    # Pool all reads from all marker files
    all_records = []
    for marker in markers:
        fq = barcode_dir / f"filtered_reads_{marker}.fastq.gz"
        if fq.exists():
            recs = read_fastq_gz(fq)
            all_records.extend(recs)

    if not all_records:
        return {"total": 0}

    # Classify each read
    bins: dict = {m: [] for m in markers}
    discarded = []
    for rec in all_records:
        seq = rec[1]
        hit = classify_read(seq, primers, max_err_rate=max_err_rate)
        if hit and hit in bins:
            bins[hit].append(rec)
        else:
            discarded.append(rec)

    # Write back (create pretrim backup first)
    for marker in markers:
        fq = barcode_dir / f"filtered_reads_{marker}.fastq.gz"
        pretrim = barcode_dir / f"filtered_reads_{marker}.pretrim.fastq.gz"
        # Backup original (pre-reclassification) before overwriting
        if fq.exists() and (not pretrim.exists() or
                             fq.stat().st_mtime > pretrim.stat().st_mtime):
            shutil.copy2(fq, pretrim)
        write_fastq_gz(bins[marker], fq)

    stats = {m: len(v) for m, v in bins.items()}
    stats["_discarded"] = len(discarded)
    stats["_total"] = len(all_records)
    return stats


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


def trim_barcode_marker(barcode_dir: Path, marker: str, fwd_primer: str,
                        rev_primer: str, cutadapt_path: str, threads: int,
                        discard_untrimmed: bool) -> tuple:
    """
    Trim one barcode/marker file with cutadapt.
    Input file is already reclassified and backed up as .pretrim.fastq.gz.
    Returns (ok: bool, message: str).
    """
    fq_gz = barcode_dir / f"filtered_reads_{marker}.fastq.gz"
    if not fq_gz.exists():
        return False, "no input file"

    n_before = count_reads(fq_gz)
    if n_before == 0:
        return True, "0 reads (no primer hits for this marker in this barcode)"

    trimmed = barcode_dir / f"filtered_reads_{marker}.trimmed.fastq.gz"

    cmd = [
        cutadapt_path,
        "-g", fwd_primer,
        "-a", rc(rev_primer),
        "--revcomp",
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
        description="Primer-based reclassification + trimming of ONT reads.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing barcode subdirectories")
    parser.add_argument("--markers", default="18S,COI",
                        help="Comma-separated marker names (default: 18S,COI)")
    parser.add_argument("--threads", type=int, default=4,
                        help="CPU cores for cutadapt (default: 4)")
    parser.add_argument("--no_reclassify", action="store_true",
                        help="Skip primer-based reclassification; only trim")
    parser.add_argument("--discard_untrimmed", action="store_true",
                        help="Discard reads where no primer detected during trimming")
    parser.add_argument("--skip", action="store_true",
                        help="Skip all processing (pass-through)")
    args = parser.parse_args()

    if args.skip:
        print("Primer trimming/reclassification skipped (--skip flag).")
        return 0

    cutadapt_path = find_tool("cutadapt") or find_tool("./env/bin/cutadapt")
    if not cutadapt_path:
        print("WARNING: cutadapt not found — skipping.", file=sys.stderr)
        print("Install with:  env/bin/pip install cutadapt", file=sys.stderr)
        return 0

    print(f"cutadapt: {cutadapt_path}")

    primers_tsv = Path(__file__).resolve().parent / "primers.tsv"
    primers     = load_primers(primers_tsv)
    print(f"Primers: {primers_tsv if primers_tsv.exists() else 'built-in defaults'}")

    markers   = [m.strip() for m in args.markers.split(",")]
    input_dir = Path(args.input_dir)

    barcode_dirs = sorted(
        d for d in input_dir.iterdir()
        if d.is_dir() and d.name.startswith("barcode")
    )
    if not barcode_dirs:
        print(f"No barcode directories found in {input_dir}")
        return 0

    # ── Phase 1: Primer-based reclassification ──────────────────────────────
    if not args.no_reclassify:
        print(f"\n=== Phase 1: Primer-based reclassification ({len(barcode_dirs)} barcodes) ===")
        print(f"  Markers: {', '.join(markers)}")
        total_in = total_kept = total_disc = 0
        for bdir in barcode_dirs:
            stats = reclassify_barcode(bdir, markers, primers)
            total_in  += stats.get("_total", 0)
            total_disc += stats.get("_discarded", 0)
            per_m = "  ".join(f"{m}:{stats.get(m, 0)}" for m in markers)
            disc  = stats.get("_discarded", 0)
            total = stats.get("_total", 0)
            total_kept += total - disc
            pct_disc = 100 * disc / max(total, 1)
            print(f"  {bdir.name}: {per_m}  discarded:{disc} ({pct_disc:.1f}%)")

        pct_total_disc = 100 * total_disc / max(total_in, 1)
        print(f"\n  Total: {total_in:,} reads in → {total_kept:,} kept, "
              f"{total_disc:,} discarded ({pct_total_disc:.1f}%)")
    else:
        print("\n=== Phase 1: Reclassification skipped (--no_reclassify) ===")

    # ── Phase 2: Cutadapt trimming ───────────────────────────────────────────
    print(f"\n=== Phase 2: Cutadapt trimming ===")
    for marker in markers:
        print(f"\n--- Marker: {marker} ---")
        if marker not in primers:
            print(f"  [SKIP] No primer sequences defined for '{marker}'.")
            continue

        fwd = primers[marker]["fwd"]
        rev = primers[marker]["rev"]
        print(f"  Fwd: {fwd}")
        print(f"  Rev: {rev}  →  RC(rev): {rc(rev)}")

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
