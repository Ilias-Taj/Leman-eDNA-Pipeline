#!/usr/bin/env python3
"""
run_validation.py

Consensus-based validation.
Updated to use 'samtools consensus' for frequency-based consensus generation.
"""

import argparse
import gzip
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
import re

def run_cmd(cmd, shell=False):
    print("Running:", cmd if shell else " ".join(map(str, cmd)))
    try:
        if shell:
            subprocess.run(cmd, check=True, shell=True, executable='/bin/bash')
        else:
            subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: command failed (exit {e.returncode}): {e}", file=sys.stderr)
        raise

def load_fasta(db_path):
    seqs = {}
    cur_id = None
    cur_seq_parts = []
    with open(db_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_seq_parts)
                header = line[1:]
                cur_id = header.split()[0]
                cur_seq_parts = []
            else:
                cur_seq_parts.append(line)
    if cur_id is not None:
        seqs[cur_id] = "".join(cur_seq_parts)
    return seqs

def parse_paf_group_reads(paf_file, mapq_threshold=40):
    groups = defaultdict(set)
    with open(paf_file, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            qname = cols[0]
            tname = cols[5]
            try:
                mapq = int(cols[11])
            except ValueError:
                mapq = 0
            if mapq >= mapq_threshold:
                groups[tname].add(qname)
    return groups

def extract_reads_by_names(reads_gz, names_set, out_fastq):
    names = set(names_set)
    found = 0
    opener = gzip.open if str(reads_gz).endswith('.gz') else open
    
    # Check file type (FASTA vs FASTQ)
    with opener(reads_gz, 'rt') as inf:
        first = inf.read(1)
    
    # Simple FASTQ extraction logic
    if first == '@':
        with opener(reads_gz, 'rt') as inf, open(out_fastq, 'wt') as outf:
            while True:
                hdr = inf.readline()
                if not hdr: break
                seq = inf.readline()
                plus = inf.readline()
                qual = inf.readline()
                if not qual: break
                
                token = hdr.strip().split()[0][1:] # Remove @
                if token in names:
                    outf.write(hdr)
                    outf.write(seq)
                    outf.write(plus)
                    outf.write(qual)
                    found += 1
                    names.remove(token)
                if not names: break
    return found

def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)

def make_single_ref_fasta(ref_id, seq, out_path):
    with open(out_path, "wt") as fh:
        fh.write(f">{ref_id}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")

def build_consensus(ref_fa, reads_fastq, work_dir, threads=4):
    """
    Uses 'samtools consensus' for a frequency-based majority vote.
    """
    ref_fa = str(ref_fa)
    reads_fastq = str(reads_fastq)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    sam = work_dir / "aln.sam"
    bam = work_dir / "aln.bam"
    sorted_bam = work_dir / "aln.sorted.bam"
    consensus = work_dir / "consensus.fa"

    # 1) minimap2 -> SAM
    cmd = ["minimap2", "-a", "-x", "map-ont", "-t", str(threads), ref_fa, reads_fastq]
    with open(sam, "wb") as samfh:
        print("Running:", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, stdout=samfh)
        except subprocess.CalledProcessError:
            print("ERROR: minimap2 failed", file=sys.stderr)
            return None

    # 2) samtools view -> bam
    run_cmd(["samtools", "view", "-bS", str(sam), "-o", str(bam)])
    
    # 3) samtools sort
    run_cmd(["samtools", "sort", "-o", str(sorted_bam), str(bam)])
    
    # 4) samtools index
    run_cmd(["samtools", "index", str(sorted_bam)])

    # 5) samtools consensus (The New Logic)
    cmd_cons = [
        "samtools", "consensus",
        "-f", "fasta",
        "-a",
        "--min-depth", "1",
        "-o", str(consensus),
        str(sorted_bam)
    ]
    
    print(f"Generating frequency-based consensus for {sorted_bam}...")
    run_cmd(cmd_cons)

    return consensus

def parse_paf_best_hit(paf_file):
    best = None
    with open(paf_file, "rt") as fh:
        for line in fh:
            if not line.strip(): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11: continue
            tname = cols[5]
            try: nmatch = int(cols[9])
            except: nmatch = 0
            try: aln_len = int(cols[10])
            except: aln_len = 0
            
            if best is None or nmatch > best[1]:
                best = (tname, nmatch, aln_len)
    return best

def main():
    parser = argparse.ArgumentParser(description="Consensus-based validation")
    parser.add_argument("--paf_file", required=True)
    parser.add_argument("--reads_file", required=True)
    parser.add_argument("--db_file", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--mapq", type=int, default=40)
    parser.add_argument("--min_reads", type=int, default=20)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    paf = Path(args.paf_file)
    reads = Path(args.reads_file)
    db = Path(args.db_file)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading reference DB...")
    refs = load_fasta(str(db))
    print("Grouping reads...")
    groups = parse_paf_group_reads(str(paf), mapq_threshold=args.mapq)

    report_rows = []
    for ref_id, read_names in groups.items():
        read_count = len(read_names)
        print(f"Reference {ref_id}: {read_count} high-MAPQ reads")
        if read_count < args.min_reads:
            report_rows.append((ref_id, read_count, 0, "", 0.0, "TOO_FEW_READS"))
            continue

        safe = safe_name(ref_id)
        reads_fastq = outdir / f"{safe}.reads.fastq"
        found = extract_reads_by_names(str(reads), read_names, str(reads_fastq))

        if ref_id not in refs:
            report_rows.append((ref_id, found, 0, "", 0.0, "NO_REF"))
            continue
        
        single_ref_fa = outdir / f"{safe}.ref.fa"
        make_single_ref_fasta(ref_id, refs[ref_id], single_ref_fa)

        work = outdir / f"work_{safe}"
        consensus_fa = build_consensus(single_ref_fa, reads_fastq, work, threads=args.threads)
        
        if consensus_fa is None or not Path(consensus_fa).exists():
            report_rows.append((ref_id, found, 0, "", 0.0, "CONS_FAIL"))
            continue

        with open(consensus_fa, "rt") as fh:
            seq = "".join([l.strip() for l in fh if not l.startswith(">")])
        cons_len = len(seq)

        paf_out = work / "consensus_vs_db.paf"
        cmd = ["minimap2", "-x", "asm5", "-t", str(args.threads), str(db), str(consensus_fa)]
        try:
            with open(paf_out, "wb") as pf:
                subprocess.run(cmd, check=True, stdout=pf)
        except:
            report_rows.append((ref_id, found, cons_len, "", 0.0, "CONS_ALIGN_FAIL"))
            continue

        best = parse_paf_best_hit(str(paf_out))
        if best is None:
            report_rows.append((ref_id, found, cons_len, "", 0.0, "NO_ALIGN"))
            continue
            
        best_tname, nmatch, aln_len = best
        identity = (nmatch / aln_len * 100.0) if aln_len > 0 else 0.0
        status = "PASS" if best_tname == ref_id else "FAIL"
        report_rows.append((ref_id, found, cons_len, best_tname, identity, status))

    report_file = outdir / "validation_report.csv"
    with open(report_file, "wt") as rf:
        rf.write("target_species,read_count,consensus_length,consensus_best_hit,consensus_identity_pct,validation_status\n")
        for row in report_rows:
            rf.write(",".join([str(x) for x in row]) + "\n")

    print(f"Wrote validation report to: {report_file}")

if __name__ == "__main__":
    main()