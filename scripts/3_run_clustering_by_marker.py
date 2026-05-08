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
        "--threads",  str(threads),
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
    """Generate per-OTU consensus sequences from isONclust3 cluster FASTQ files.

    For each cluster:
      - Converts cluster FASTQ to FASTA
      - Runs vsearch --cluster_fast at 80 % identity to produce a consensus
      - Writes >centroid={OTU_id};seqs={n_reads} to the output FASTA

    Args:
        isonclust_out_dir: Path to the clustering/ directory returned by
                           run_isonclust3 (contains final_clusters.tsv and
                           fastq_files/).
        working_dir:       temp_clustering/ directory for output files.
        marker:            Marker name (18S, COI, JEDI).
        vsearch_path:      Path to vsearch executable.
        min_cluster_size:  Skip clusters with fewer reads than this.

    Returns:
        Tuple (consensus_fasta_path, cluster_to_otu_dict) where
        cluster_to_otu_dict maps isONclust3 cluster_id -> OTU ID string.
    """
    isonclust_out_dir = Path(isonclust_out_dir)
    working_dir       = Path(working_dir)

    clusters_tsv = isonclust_out_dir / "final_clusters.tsv"
    if not clusters_tsv.exists():
        print(f"ERROR: {clusters_tsv} not found", file=sys.stderr)
        sys.exit(1)

    # Parse cluster_id -> [read_ids]
    cluster_reads = {}
    with open(clusters_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            cid, rid = parts[0], parts[1]
            cluster_reads.setdefault(cid, []).append(rid)

    def _sort_key(x):
        try:
            return int(x)
        except ValueError:
            return x

    sorted_ids = sorted(cluster_reads.keys(), key=_sort_key)
    fastq_dir  = isonclust_out_dir / "fastq_files"
    consensus_fasta = working_dir / f"consensus_{marker}.fasta"

    cluster_to_otu = {}
    otu_counter    = 0
    total_seqs     = 0

    with open(consensus_fasta, "w") as out_f:
        for cluster_id in sorted_ids:
            read_ids = cluster_reads[cluster_id]
            n_reads  = len(read_ids)

            if n_reads < min_cluster_size:
                continue

            cluster_fastq = fastq_dir / f"{cluster_id}.fastq"
            if not cluster_fastq.exists():
                continue

            otu_counter += 1
            otu_id = f"OTU_{marker}_{otu_counter:06d}"
            cluster_to_otu[cluster_id] = otu_id

            # Convert cluster FASTQ -> FASTA (simple 4-line parser)
            tmp_fasta = working_dir / f"_tmp_cluster_{cluster_id}.fasta"
            with open(tmp_fasta, "w") as fa_out:
                with open(cluster_fastq, "r") as fq_in:
                    buf = []
                    for line in fq_in:
                        buf.append(line.rstrip("\n"))
                        if len(buf) == 4:
                            fa_out.write(f">{buf[0][1:]}\n{buf[1]}\n")
                            buf = []

            if n_reads == 1:
                # Single-read cluster (shouldn't occur with min_cluster_size=2)
                with open(tmp_fasta, "r") as f_in:
                    lines_list = [l.rstrip("\n") for l in f_in]
                if len(lines_list) >= 2:
                    out_f.write(
                        f">centroid={otu_id};seqs={n_reads}\n"
                        f"{lines_list[1]}\n"
                    )
                    total_seqs += 1
                tmp_fasta.unlink(missing_ok=True)
                continue

            # Run vsearch to get a consensus of reads within the cluster
            tmp_cons = working_dir / f"_tmp_cons_{cluster_id}.fasta"
            cmd = [
                vsearch_path,
                "--cluster_fast", str(tmp_fasta),
                "--id",           "0.80",
                "--consout",      str(tmp_cons),
                "--minseqlength", "100",
                "--quiet",
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            tmp_fasta.unlink(missing_ok=True)

            if res.returncode != 0 or not tmp_cons.exists():
                continue

            # Take only the first consensus sequence
            with open(tmp_cons, "r") as cons_in:
                found_seq = False
                for line in cons_in:
                    if line.startswith(">"):
                        found_seq = True
                        continue
                    if found_seq and line.strip():
                        out_f.write(
                            f">centroid={otu_id};seqs={n_reads}\n"
                            f"{line.strip()}\n"
                        )
                        total_seqs += 1
                        break

            tmp_cons.unlink(missing_ok=True)

    print(
        f"Generated {total_seqs} {marker} consensus sequences "
        f"from {len(cluster_reads)} isONclust3 clusters"
    )
    return consensus_fasta, cluster_to_otu


# ── Assignment generation ─────────────────────────────────────────────────────

def parse_isonclust_to_assignment(isonclust_out_dir, output_file, marker,
                                   cluster_to_otu, chimera_ids=None):
    """Write a per-read OTU assignment file from isONclust3 clusters.

    Reads final_clusters.tsv and maps each read to its OTU using the
    cluster_to_otu dict produced by generate_consensus_from_clusters.
    Chimeric OTUs are excluded.

    Output format (one line per read):
        {uuid}  {OTU_id}  {barcode}  {marker}

    Args:
        isonclust_out_dir: clustering/ directory from run_isonclust3.
        output_file:       Destination .txt file.
        marker:            Marker name.
        cluster_to_otu:    dict {cluster_id -> OTU_id} (from
                           generate_consensus_from_clusters).
        chimera_ids:       Set of OTU IDs identified as chimeric.
    """
    isonclust_out_dir = Path(isonclust_out_dir)
    chimera_ids       = chimera_ids or set()
    clusters_tsv      = isonclust_out_dir / "final_clusters.tsv"

    written_count = 0
    with open(output_file, "w") as f_out:
        with open(clusters_tsv, "r") as f_in:
            for line in f_in:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                cluster_id, read_full_id = parts[0], parts[1]

                otu_id = cluster_to_otu.get(cluster_id)
                if otu_id is None:
                    continue  # cluster filtered (too small or missing fastq)
                if otu_id in chimera_ids:
                    continue  # chimeric OTU

                # Parse read ID: {uuid}|{barcode}|{marker}
                id_parts = read_full_id.split("|")
                if len(id_parts) < 3:
                    continue
                uuid    = id_parts[0]
                barcode = id_parts[1]
                m       = id_parts[2]

                f_out.write(f"{uuid}\t{otu_id}\t{barcode}\t{m}\n")
                written_count += 1

    n_otus = len(set(cluster_to_otu.values()))
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