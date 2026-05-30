#!/usr/bin/env python3
"""
2_classify_markers.py - CLASSIFY READS BY MARKER

Separates reads into marker bins using two strategies:
1. minimap2-based classification (primary, if --use_minimap2 is set)
2. Length-based fallback (always active for unclassified reads)

Length windows:
- 18S rRNA: 1,500-2,800 bp
- COI:      500-900 bp   (standard Folmer/Leray primers, ~658 bp amplicon)
- JEDI:     250-500 bp   (JEDI universal rRNA V4-V5 primers, 515F-Y/926R)

Use --markers to select which markers to search for (default: 18S,COI).
Reads outside all selected marker ranges are discarded as noise.

Prerequisite: Run 1_run_preprocessing.py first

For each barcode's filtered_reads.fastq.gz, outputs:
- filtered_reads_{marker}.fastq.gz  (one per selected marker)
- marker_counts.txt

Usage:
    # Water data (18S + COI) with minimap2 classification
    python scripts/2_classify_markers.py \
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \
        --markers 18S,COI \
        --use_minimap2 --db_18S pr2 --db_COI midori2

    # Soil data (JEDI + COI) length-only
    python scripts/2_classify_markers.py \
        --input_dir out/Soil_eDNA_JEDI_COI_14_01_26 --markers JEDI,COI

Changelog:
    2025-05-12  Hardcode minimap2 threads to 4 (sufficient for length-based marker classification).
"""

import argparse
import gzip
import subprocess
import sys
from pathlib import Path

# Ensure scripts/ is on the import path so utils.py can be found from any cwd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import SKIP_DIRS, find_tool, find_minimap2

# Sequence length thresholds (bp) for each marker
MARKER_RANGES = {
    "18S":  (1500, 2800),   # 18S rRNA full-length (~1.8 kb amplicon)
    "COI":  (500,  900),    # COI Folmer/Leray (~658 bp target)
    "JEDI": (250,  500),    # JEDI rRNA V4-V5, 515F-Y/926R (~390 bp target)
}

VALID_MARKERS = list(MARKER_RANGES.keys())

# Short name → filename inside refs/
_DB_NAME_MAP = {
    "pr2":     "pr2_18S.udb",
    "silva":   "silva_18S.udb",
    "midori2": "midori2_COI.udb",
    "midori":  "midori2_COI.udb",
    "ekoi":    "eKOI_COI.udb",
    "porter":  "porter_COI.udb",
}

# Auto-detect fallback order when no DB is specified
_MARKER_DEFAULT_DBS = {
    "18S":  ["pr2_18S.udb", "silva_18S.udb",
             "pr2_version_5.1.1_SSU_dada2.fasta.gz"],
    "COI":  ["midori2_COI.udb", "eKOI_COI.udb", "porter_COI.udb",
             "eKOI_COI_SINTAX.fasta"],
    "JEDI": ["pr2_18S.udb", "silva_18S.udb"],
}


# ── Database helpers ─────────────────────────────────────────────────────────

def resolve_db_path(db_name_or_path, refs_dir, marker):
    """Return a Path to the reference database, or None if not found."""
    refs_dir = Path(refs_dir)
    if db_name_or_path:
        p = Path(db_name_or_path)
        if p.exists():
            return p
        key = db_name_or_path.lower()
        if key in _DB_NAME_MAP:
            candidate = refs_dir / _DB_NAME_MAP[key]
            if candidate.exists():
                return candidate
    # Auto-detect from defaults
    for fname in _MARKER_DEFAULT_DBS.get(marker, []):
        p = refs_dir / fname
        if p.exists():
            return p
    return None


def _read_fasta_sequences(path):
    """Read (header, seq) pairs from a plain FASTA file."""
    sequences = []
    cur_hdr, cur_seq = None, []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_hdr is not None:
                    sequences.append((cur_hdr, "".join(cur_seq)))
                cur_hdr = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
    if cur_hdr is not None:
        sequences.append((cur_hdr, "".join(cur_seq)))
    return sequences


def _read_fasta_sequences_gz(path):
    """Read (header, seq) pairs from a gzipped FASTA file."""
    sequences = []
    cur_hdr, cur_seq = None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_hdr is not None:
                    sequences.append((cur_hdr, "".join(cur_seq)))
                cur_hdr = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
    if cur_hdr is not None:
        sequences.append((cur_hdr, "".join(cur_seq)))
    return sequences


# ── Reference extraction ─────────────────────────────────────────────────────

def extract_marker_references(refs_dir, marker, db_path_or_name, vsearch_path,
                               n_refs=50, length_range=None):
    """Extract representative reference sequences for use with minimap2.

    Args:
        refs_dir:        Path to refs/ directory.
        marker:          Marker name (18S, COI, JEDI).
        db_path_or_name: Short name (pr2, silva, …) or path to .udb/.fasta[.gz].
        vsearch_path:    Path to vsearch executable.
        n_refs:          Max number of reference sequences to extract.
        length_range:    (lo, hi) bp; defaults to MARKER_RANGES[marker].

    Returns:
        Path to refs/minimap2_refs/{marker}.fasta, or None on failure.
    """
    refs_dir = Path(refs_dir)
    out_dir = refs_dir / "minimap2_refs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fasta = out_dir / f"{marker}.fasta"

    # Use cached ref file if it already exists (DB rarely changes between runs)
    if out_fasta.exists() and out_fasta.stat().st_size > 0:
        print(f"  [REF] Using cached {marker} references: {out_fasta}")
        return out_fasta

    if length_range is None:
        length_range = MARKER_RANGES.get(marker, (0, 99999))
    lo, hi = length_range

    db_path = resolve_db_path(db_path_or_name, refs_dir, marker)
    if db_path is None:
        print(f"  [WARN] No reference database found for {marker}", file=sys.stderr)
        return None

    suffix = db_path.suffix.lower()
    name_lower = db_path.name.lower()

    sequences = []

    if suffix == ".udb":
        # Stream udb2fasta output through a pipe to avoid writing hundreds of
        # thousands of sequences to disk before length filtering.
        import subprocess as _sp
        proc = _sp.Popen(
            [vsearch_path, "--udb2fasta", str(db_path),
             "--output", "-", "--quiet"],
            stdout=_sp.PIPE, stderr=_sp.PIPE, text=True,
        )
        sequences = []
        hdr, seq = "", ""
        for line in proc.stdout:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if hdr and lo <= len(seq) <= hi:
                    sequences.append((hdr, seq))
                    if len(sequences) >= n_refs * 20:  # keep a generous pool
                        break
                hdr, seq = line[1:], ""
            else:
                seq += line
        if hdr and lo <= len(seq) <= hi:
            sequences.append((hdr, seq))
        proc.stdout.close()
        proc.wait()
        if proc.returncode not in (0, None) and not sequences:
            err = proc.stderr.read()[:200] if proc.stderr else ""
            print(f"  [WARN] vsearch udb2fasta failed for {marker}: {err}",
                  file=sys.stderr)
            return None

    elif name_lower.endswith(".fasta.gz") or name_lower.endswith(".fa.gz"):
        sequences = _read_fasta_sequences_gz(db_path)

    elif suffix in (".fasta", ".fa", ".fas"):
        sequences = _read_fasta_sequences(db_path)

    else:
        print(f"  [WARN] Unsupported database format: {db_path}", file=sys.stderr)
        return None

    filtered = [(h, s) for h, s in sequences if lo <= len(s) <= hi]
    if not filtered:
        print(f"  [WARN] No sequences in length range {lo}-{hi} bp "
              f"found in {db_path.name}", file=sys.stderr)
        return None

    # Sample evenly up to n_refs
    if len(filtered) > n_refs:
        step = max(1, len(filtered) // n_refs)
        filtered = filtered[::step][:n_refs]

    # Write with marker-prefixed headers so PAF parsing can identify the marker
    with open(out_fasta, "w") as f:
        for header, seq in filtered:
            safe = header.split()[0].replace("|", "_").replace(";", "_")
            f.write(f">{marker}_{safe}\n{seq}\n")

    print(f"  [REF] Extracted {len(filtered)} {marker} references "
          f"from {db_path.name}")
    return out_fasta


# ── minimap2 classification ───────────────────────────────────────────────────

def classify_reads_minimap2(fastq_gz, refs_per_marker, minimap2_path,
                             temp_dir, active_markers, threads=4):
    """Classify reads via minimap2 against combined marker reference sequences.

    Reference headers must start with '{marker}_' (written by
    extract_marker_references) so the marker can be recovered from the
    target name in the PAF output.

    Args:
        fastq_gz:        Path to the barcode's filtered_reads.fastq.gz.
        refs_per_marker: dict {marker: Path} from extract_marker_references.
        minimap2_path:   Path to minimap2 executable.
        temp_dir:        Directory for temporary PAF / combined-ref files.
        active_markers:  List of marker names to consider.

    Returns:
        dict mapping read_id (str) -> marker (str) for confidently mapped reads.
    """
    temp_dir = Path(temp_dir)

    # Build combined reference FASTA
    combined_refs = temp_dir / "minimap2_combined_refs.fasta"
    has_refs = False
    with open(combined_refs, "w") as out_f:
        for marker in active_markers:
            ref_path = refs_per_marker.get(marker)
            if ref_path and Path(ref_path).exists():
                with open(ref_path, "r") as in_f:
                    out_f.write(in_f.read())
                has_refs = True

    if not has_refs:
        combined_refs.unlink(missing_ok=True)
        return {}

    paf_file = temp_dir / "minimap2_classification.paf"
    cmd = [minimap2_path, "-x", "map-ont", "--secondary=no",
           "--no-long-join", "-t", str(threads),
           str(combined_refs), str(fastq_gz)]

    with open(paf_file, "w") as paf_out:
        result = subprocess.run(cmd, stdout=paf_out, stderr=subprocess.PIPE)
    combined_refs.unlink(missing_ok=True)

    if result.returncode != 0:
        print(f"  [WARN] minimap2 failed: {result.stderr.decode()[:300]}",
              file=sys.stderr)
        paf_file.unlink(missing_ok=True)
        return {}

    # PAF columns: qname qlen qstart qend strand tname tlen tstart tend
    #              nmatch alen mapq ...
    read_to_marker = {}
    with open(paf_file, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue
            read_id  = cols[0]
            read_len = int(cols[1])
            aln_len  = int(cols[10])
            mapq     = int(cols[11])
            ref_name = cols[5]

            if mapq > 20 and aln_len > 0.7 * read_len:
                marker = ref_name.split("_")[0]
                if marker in active_markers:
                    read_to_marker[read_id] = marker

    paf_file.unlink(missing_ok=True)
    return read_to_marker


# ── Length-based fallback ─────────────────────────────────────────────────────

def classify_read(seq, read_id, active_markers):
    """Classify a single read by sequence length. Returns marker or 'ambiguous'."""
    seq_len = len(seq)
    for marker in active_markers:
        lo, hi = MARKER_RANGES[marker]
        if lo <= seq_len <= hi:
            return marker
    return "ambiguous"


# ── Per-barcode processing ────────────────────────────────────────────────────

def process_barcode(barcode_dir, active_markers,
                    use_minimap2=False, refs_per_marker=None, minimap2_path=None,
                    threads=4):
    """Split filtered_reads.fastq.gz into per-marker FASTQ files.

    Uses minimap2 as the primary classifier (when requested) and falls back
    to length-based classification for unclassified reads.
    """
    barcode_name = barcode_dir.name
    filtered_reads = barcode_dir / "filtered_reads.fastq.gz"

    if not filtered_reads.exists():
        print(f"  [SKIP] {barcode_name} - no filtered_reads.fastq.gz",
              file=sys.stderr)
        return None

    # Run minimap2 pre-classification if requested
    minimap2_classifications = {}
    if use_minimap2 and refs_per_marker and minimap2_path:
        try:
            minimap2_classifications = classify_reads_minimap2(
                filtered_reads, refs_per_marker, minimap2_path,
                barcode_dir, active_markers, threads=threads,
            )
            if minimap2_classifications:
                print(f"  [minimap2] {barcode_name}: "
                      f"{len(minimap2_classifications)} reads pre-classified")
        except Exception as e:
            print(f"  [WARN] minimap2 failed for {barcode_name}: {e}",
                  file=sys.stderr)
            minimap2_classifications = {}

    # Open per-marker output handles
    output_handles = {}
    for marker in active_markers:
        output_handles[marker] = gzip.open(
            barcode_dir / f"filtered_reads_{marker}.fastq.gz", "wt"
        )

    counts = {m: 0 for m in active_markers}
    counts["ambiguous"] = 0

    with gzip.open(filtered_reads, "rt") as in_f:
        lines = []
        for line in in_f:
            lines.append(line)
            if len(lines) == 4:
                header, seq, plus, qual = lines
                read_id = header[1:].strip().split()[0]

                # Primary: minimap2; fallback: length-based
                marker = minimap2_classifications.get(read_id)
                if marker is None or marker not in active_markers:
                    marker = classify_read(seq.strip(), read_id, active_markers)

                if marker in output_handles:
                    output_handles[marker].write(header)
                    output_handles[marker].write(seq)
                    output_handles[marker].write(plus)
                    output_handles[marker].write(qual)
                    counts[marker] += 1
                else:
                    counts["ambiguous"] += 1

                lines = []

    for fh in output_handles.values():
        fh.close()

    total = sum(counts.values())
    parts = ", ".join(f"{m}={counts[m]}" for m in active_markers)
    print(f"{barcode_name}: {parts}, ambiguous={counts['ambiguous']}, total={total}")

    counts_file = barcode_dir / "marker_counts.txt"
    with open(counts_file, "w") as f:
        for m in active_markers:
            f.write(f"{m}: {counts[m]}\n")
        f.write(f"Ambiguous: {counts['ambiguous']}\n")
        f.write(f"Total: {total}\n")

    return {
        "barcode": barcode_name,
        **{m: counts[m] for m in active_markers},
        "ambiguous": counts["ambiguous"],
        "total": total,
    }


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Classify eDNA reads by marker (minimap2 + length fallback)"
    )
    parser.add_argument(
        "--input_dir", required=True,
        help="Input directory containing barcode folders "
             "(e.g., out/Water_eDNA_18S_COI_14_01_26)"
    )
    parser.add_argument(
        "--markers", default="18S,COI",
        help=(
            "Comma-separated list of markers to classify. "
            f"Valid markers: {', '.join(VALID_MARKERS)}. "
            "Default: 18S,COI. For soil JEDI data use: JEDI,COI"
        )
    )
    parser.add_argument(
        "--use_minimap2", action="store_true", default=False,
        help="Use minimap2 as primary classifier (falls back to length-based)"
    )
    parser.add_argument(
        "--refs_dir", default="refs/",
        help="Directory containing reference databases (default: refs/)"
    )
    parser.add_argument(
        "--db_18S", default="",
        help="18S reference DB: pr2, silva, or path to .udb/.fasta[.gz]"
    )
    parser.add_argument(
        "--db_COI", default="",
        help="COI reference DB: midori2, ekoi, porter, or path"
    )
    parser.add_argument(
        "--db_JEDI", default="",
        help="JEDI reference DB: pr2, silva, or path (JEDI uses 18S databases)"
    )

    args = parser.parse_args()

    # Validate markers
    active_markers = [m.strip().upper() for m in args.markers.split(",")]
    for m in active_markers:
        if m not in VALID_MARKERS:
            print(f"ERROR: Unknown marker '{m}'. "
                  f"Valid markers: {', '.join(VALID_MARKERS)}", file=sys.stderr)
            sys.exit(1)

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"ERROR: Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    refs_dir = Path(args.refs_dir)

    print(f"Classifying reads by marker ({', '.join(active_markers)})...")
    print(f"Input directory: {input_dir}")
    print("Length ranges:")
    for m in active_markers:
        lo, hi = MARKER_RANGES[m]
        print(f"  {m}: {lo}-{hi} bp")
    print()

    # Prepare minimap2 resources if requested
    refs_per_marker = {}
    minimap2_path = None
    if args.use_minimap2:
        minimap2_path = find_minimap2()
        vsearch_path = find_tool("vsearch")
        db_args = {"18S": args.db_18S, "COI": args.db_COI, "JEDI": args.db_JEDI}
        print("Extracting minimap2 reference sequences...")
        for marker in active_markers:
            db = db_args.get(marker) or None
            ref_path = extract_marker_references(
                refs_dir, marker, db, vsearch_path, n_refs=50,
            )
            if ref_path:
                refs_per_marker[marker] = ref_path
        if not refs_per_marker:
            print("[WARN] No minimap2 references extracted; "
                  "using length-only classification", file=sys.stderr)
        print()

    # Process each barcode directory
    results = []
    for barcode_dir in sorted(input_dir.iterdir()):
        if not barcode_dir.is_dir():
            continue
        if barcode_dir.name in SKIP_DIRS:
            continue

        result = process_barcode(
            barcode_dir, active_markers,
            use_minimap2=args.use_minimap2,
            refs_per_marker=refs_per_marker,
            minimap2_path=minimap2_path,
            threads=4,
        )
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