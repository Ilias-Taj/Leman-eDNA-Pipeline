#!/usr/bin/env python3
"""
3_run_clustering_by_marker.py - CLUSTERING + CHIMERA DETECTION

Clusters reads per marker using isONclust3, generates per-cluster consensus
sequences with VSEARCH, then removes chimeric OTUs.

Prerequisite: Run 2_classify_markers.py first

Usage:
    python scripts/3_run_clustering_by_marker.py \
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \
        --output_dir out/Water_eDNA_18S_COI_14_01_26 \
        --threads 14 \
        --isonclust3 ./tools/isONclust3/target/release/isONclust3

Output:
    Creates temp_clustering/ folder with:
    - consensus_{marker}_clean.fasta  (chimera-filtered consensus sequences)
    - global_otu_assignment.txt       (read-to-OTU mapping across all markers)
    - global_otu_assignment_{marker}.txt  (per-marker assignment files)
"""

import argparse
import gzip
import subprocess
import sys
import time
from pathlib import Path

# Ensure scripts/ is on the import path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import find_tool, find_isonclust3, timer, SKIP_DIRS


# ── Chimera detection (UNCHANGED) ────────────────────────────────────────────

def remove_chimeras(working_dir, marker, vsearch_path):
    """Detect and remove chimeric OTUs using VSEARCH --uchime_denovo."""
    consensus_in    = working_dir / f"consensus_{marker}.fasta"
    consensus_clean = working_dir / f"consensus_{marker}_clean.fasta"
    chimeras_file   = working_dir / f"chimeras_{marker}.fasta"

    cmd = [
        vsearch_path,
        "--uchime_denovo", str(consensus_in),
        "--nonchimeras",   str(consensus_clean),
        "--chimeras",      str(chimeras_file),
        "--xsize",
    ]

    print(f"Running Chimera Detection for {marker}...")
    subprocess.run(cmd, check=True)
    return consensus_clean, chimeras_file


def load_chimera_ids(chimeras_fasta):
    """Load chimera centroid IDs from a FASTA file.

    VSEARCH consensus headers look like:
        >centroid=OTU_18S_000001;seqs=N
    We strip 'centroid=' so the IDs match OTU ID format used in assignments.
    """
    chimera_ids = set()
    if not chimeras_fasta.exists():
        return chimera_ids
    with open(chimeras_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                header   = line[1:].strip().split()[0]
                clean_id = header.split(";")[0]
                clean_id = clean_id.replace("centroid=", "")
                chimera_ids.add(clean_id)
    return chimera_ids


# ── Read concatenation ────────────────────────────────────────────────────────

def concatenate_reads_fastq_by_marker(input_dir, marker, output_fastq):
    """Concatenate all per-barcode marker FASTQ files into one FASTQ.

    Read IDs are rewritten as  {uuid}|{barcode}|{marker}  so provenance is
    preserved through clustering and assignment.
    """
    print(f"Collecting {marker} reads from {input_dir}...")
    total_reads = 0

    with open(output_fastq, "w") as out_f:
        for sample_dir in sorted(input_dir.iterdir()):
            if not sample_dir.is_dir() or sample_dir.name in SKIP_DIRS:
                continue
            fq_path = sample_dir / f"filtered_reads_{marker}.fastq.gz"
            if not fq_path.exists():
                continue

            barcode = sample_dir.name
            with gzip.open(fq_path, "rt") as in_f:
                buf = []
                for line in in_f:
                    buf.append(line.rstrip("\n"))
                    if len(buf) == 4:
                        orig_id = buf[0][1:].split()[0]
                        seq     = buf[1]
                        qual    = buf[3]
                        out_f.write(
                            f"@{orig_id}|{barcode}|{marker}\n"
                            f"{seq}\n+\n{qual}\n"
                        )
                        total_reads += 1
                        buf = []

    print(f"Concatenated {total_reads} {marker} reads into {output_fastq}")
    return total_reads


# ── isONclust3 clustering ─────────────────────────────────────────────────────

def run_isonclust3(input_fastq, working_dir, marker, min_cluster_size,
                   threads, isonclust3_path):
    """Run isONclust3 on the combined marker FASTQ.

    Returns:
        Path to the clustering/ subdirectory inside the isONclust3 output folder.
    """
    outfolder = Path(working_dir) / f"isonclust_{marker}"
    cmd = [
        isonclust3_path,
        "--fastq",    str(input_fastq),
        "--mode",     "ont",
        "--outfolder", str(outfolder),
        "--n",        str(min_cluster_size),
        "--seeding",  "minimizer",
    ]

    print(f"Running isONclust3 for {marker}: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"isONclust3 failed for {marker}:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    clustering_dir = outfolder / "clustering"
    if not clustering_dir.exists():
        print(f"ERROR: expected isONclust3 output dir not found: {clustering_dir}",
              file=sys.stderr)
        sys.exit(1)

    return clustering_dir


# ── Consensus generation ──────────────────────────────────────────────────────

def generate_consensus_from_clusters(isonclust_out_dir, working_dir, marker,
                                     vsearch_path, min_cluster_size=2):
    """Generate per-OTU representative sequences from isONclust3 cluster FASTQ files.

    Iterates directly over the per-cluster FASTQ files in fastq_files/.
    For each cluster, picks the longest read as the OTU representative —
    longer reads are more likely to span the full amplicon.

    NOTE: The TSV cluster IDs do NOT map 1-to-1 to fastq file numbers
    (isONclust3 renumbers clusters internally), so the TSV is only used for
    the read->cluster assignment in parse_isonclust_to_assignment, not here.

    Args:
        isonclust_out_dir: Path to the clustering/ directory returned by
                           run_isonclust3 (contains fastq_files/).
        working_dir:       temp_clustering/ directory for output files.
        marker:            Marker name (18S, COI, JEDI).
        vsearch_path:      Unused — kept for API compatibility with chimera step.
        min_cluster_size:  Skip clusters with fewer reads than this.

    Returns:
        Tuple (consensus_fasta_path, read_to_otu_dict) where read_to_otu_dict
        maps each read_id to the OTU_id of its cluster.
    """
    isonclust_out_dir = Path(isonclust_out_dir)
    working_dir       = Path(working_dir)

    fastq_dir = isonclust_out_dir / "fastq_files"
    if not fastq_dir.exists():
        print(f"ERROR: fastq_files dir not found: {fastq_dir}", file=sys.stderr)
        sys.exit(1)

    # Collect all per-cluster FASTQ files and sort numerically
    fastq_files = sorted(fastq_dir.glob("*.fastq"),
                         key=lambda p: int(p.stem) if p.stem.isdigit() else p.stem)

    consensus_fasta = working_dir / f"consensus_{marker}.fasta"
    read_to_otu     = {}   # {read_id: otu_id} for all reads in kept clusters
    otu_counter     = 0
    total_seqs      = 0

    with open(consensus_fasta, "w") as out_f:
        for fq_path in fastq_files:
            # Parse all reads from this cluster FASTQ (4-line format)
            reads = []  # list of (read_id, seq) tuples
            with open(fq_path, "r") as fq_in:
                buf = []
                for line in fq_in:
                    buf.append(line.rstrip("\n"))
                    if len(buf) == 4:
                        reads.append((buf[0][1:], buf[1]))  # strip '@'
                        buf = []

            n_reads = len(reads)
            if n_reads < min_cluster_size:
                continue

            # Pick the longest read as OTU representative (most likely full-length)
            rep_read_id, rep_seq = max(reads, key=lambda r: len(r[1]))

            otu_counter += 1
            otu_id = f"OTU_{marker}_{otu_counter:06d}"

            for read_id, _ in reads:
                read_to_otu[read_id] = otu_id

            out_f.write(f">centroid={otu_id};seqs={n_reads}\n{rep_seq}\n")
            total_seqs += 1

    print(
        f"Generated {total_seqs} {marker} representative sequences "
        f"from {len(fastq_files)} isONclust3 cluster files"
    )
    return consensus_fasta, read_to_otu


# ── Assignment generation ─────────────────────────────────────────────────────

def parse_isonclust_to_assignment(isonclust_out_dir, output_file, marker,
                                   cluster_to_otu, chimera_ids=None):
    """Write a per-read OTU assignment file from isONclust3 clusters.

    Uses the read_to_otu dict (produced by generate_consensus_from_clusters)
    to map each read directly to its OTU without relying on TSV cluster IDs
    (which do not match fastq file numbers in isONclust3 output).

    Args:
        isonclust_out_dir: clustering/ directory from run_isonclust3 (unused,
                           kept for API compatibility).
        output_file:       Destination .txt file.
        marker:            Marker name.
        cluster_to_otu:    dict {read_id -> OTU_id} returned by
                           generate_consensus_from_clusters.
        chimera_ids:       Set of OTU IDs identified as chimeric.
    """
    chimera_ids = chimera_ids or set()
    read_to_otu = cluster_to_otu   # renamed arg for clarity

    written_count = 0
    with open(output_file, "w") as f_out:
        for read_full_id, otu_id in read_to_otu.items():
            if otu_id in chimera_ids:
                continue

            # Parse read ID: {uuid}|{barcode}|{marker}
            id_parts = read_full_id.split("|")
            if len(id_parts) < 3:
                continue
            uuid    = id_parts[0]
            barcode = id_parts[1]
            m       = id_parts[2]

            f_out.write(f"{uuid}\t{otu_id}\t{barcode}\t{m}\n")
            written_count += 1

    n_otus = len(set(otu_id for otu_id in read_to_otu.values()
                     if otu_id not in chimera_ids))
    print(
        f"Generated assignments for {n_otus} {marker} OTUs, "
        f"wrote {written_count} assignments."
    )
    return n_otus

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Global OTU clustering by marker using isONclust3 + VSEARCH"
    )
    parser.add_argument("--input_dir",  required=True,
                        help="Input directory containing barcode subdirectories")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for assignment files")
    parser.add_argument("--threads", type=int, default=4,
                        help="Threads for isONclust3 (default: 4)")
    parser.add_argument("--markers", default="18S,COI",
                        help="Comma-separated list of markers (default: 18S,COI)")
    parser.add_argument(
        "--isonclust3",
        default="./tools/isONclust3/target/release/isONclust3",
        help="Path to isONclust3 binary "
             "(default: ./tools/isONclust3/target/release/isONclust3)",
    )
    parser.add_argument("--min_cluster_size", type=int, default=2,
                        help="Minimum reads per cluster (default: 2)")

    args = parser.parse_args()

    active_markers = [m.strip().upper() for m in args.markers.split(",")]

    vsearch_path = find_tool("vsearch")

    # Resolve isONclust3 path
    isonclust3_path = args.isonclust3
    if not Path(isonclust3_path).exists():
        isonclust3_path = find_isonclust3()

    input_dir = Path(args.input_dir)
    out_dir   = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    temp_dir = out_dir / "temp_clustering"
    temp_dir.mkdir(exist_ok=True)

    all_otus  = {}
    timings   = {}
    total_start = time.time()

    for marker in active_markers:
        print(f"\n{'='*60}")
        print(f"PROCESSING MARKER: {marker}")
        print(f"{'='*60}")

        merged_fastq = temp_dir / f"all_reads_{marker}.fastq"

        # 1. Concatenate reads from all barcodes into one FASTQ
        with timer(f"{marker} - Concatenate reads") as t:
            num_reads = concatenate_reads_fastq_by_marker(
                input_dir, marker, merged_fastq
            )
        timings[f"{marker}_concat"] = t.elapsed

        if num_reads == 0:
            print(f"[SKIP] No {marker} reads found")
            continue

        # 2. Cluster with isONclust3
        with timer(f"{marker} - isONclust3 clustering") as t:
            isonclust_clustering_dir = run_isonclust3(
                merged_fastq, temp_dir, marker,
                args.min_cluster_size, args.threads, isonclust3_path,
            )
        timings[f"{marker}_clustering"] = t.elapsed

        # 3. Generate consensus sequences from clusters
        with timer(f"{marker} - Consensus generation") as t:
            consensus_fasta, cluster_to_otu = generate_consensus_from_clusters(
                isonclust_clustering_dir, temp_dir, marker, vsearch_path,
                args.min_cluster_size,
            )
        timings[f"{marker}_consensus"] = t.elapsed

        # 4. Chimera detection on consensus sequences
        with timer(f"{marker} - Chimera detection") as t:
            consensus_clean, chimeras_file = remove_chimeras(
                temp_dir, marker, vsearch_path
            )
            chimera_ids = load_chimera_ids(chimeras_file)
        timings[f"{marker}_chimera"] = t.elapsed

        total_consensus = sum(
            1 for l in open(consensus_fasta) if l.startswith(">")
        )
        n_chimeras = len(chimera_ids)
        n_clean    = sum(1 for l in open(consensus_clean) if l.startswith(">"))
        pct = 100 * n_chimeras / total_consensus if total_consensus else 0
        print(
            f"  Chimeras: {n_chimeras:,}/{total_consensus:,} ({pct:.1f}%) "
            f"removed, {n_clean:,} clean OTUs retained"
        )

        # 5. Build read-to-OTU assignment file (skipping chimeric OTUs)
        final_assignment = out_dir / f"global_otu_assignment_{marker}.txt"
        with timer(f"{marker} - Parse assignments") as t:
            num_otus = parse_isonclust_to_assignment(
                isonclust_clustering_dir, final_assignment, marker,
                cluster_to_otu, chimera_ids=chimera_ids,
            )
        timings[f"{marker}_parse"] = t.elapsed

        all_otus[marker] = num_otus
        print(f"[OK] {marker} OTU assignment file: {final_assignment}")

    # Merge all marker assignments into one global file
    print(f"\n{'='*60}")
    print("MERGING MARKER ASSIGNMENTS")
    print(f"{'='*60}")

    merged_assignment = out_dir / "global_otu_assignment.txt"
    with open(merged_assignment, "w") as out_f:
        out_f.write("read_name\totu_id\tbarcode\tmarker\n")
        for marker in active_markers:
            assignment_file = out_dir / f"global_otu_assignment_{marker}.txt"
            if assignment_file.exists():
                with open(assignment_file, "r") as in_f:
                    for line in in_f:
                        out_f.write(line)

    total_elapsed = time.time() - total_start

    print(f"\n{'='*60}")
    print("TIMING SUMMARY")
    print(f"{'='*60}")
    for key, elapsed in timings.items():
        print(f"  {key:30s} {elapsed:8.1f}s")
    print(f"  {'TOTAL':30s} {total_elapsed:8.1f}s")

    for marker in active_markers:
        print(f"{marker} OTUs: {all_otus.get(marker, 0)}")
    print(f"Total OTUs: {sum(all_otus.values())}")
    print(f"Global assignment file (with markers): {merged_assignment}")


if __name__ == "__main__":
    main()