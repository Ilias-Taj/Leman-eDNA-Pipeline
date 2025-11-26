#!/usr/bin/env python3
"""
5_run_validation.py

Consensus-based validation inspired by ONtoBAR.

Usage example:
  python 5_run_validation.py \
    --paf_file out/alignments.paf \
    --reads_file out/filtered_reads.fastq.gz \
    --db_file refs/sanger_db.fa \
    --output_dir out/validation

This script:
 - Parses the PAF, groups reads by reference (MAPQ filter)
 - Extracts reads for each reference with enough support
 - Builds a consensus per-reference using minimap2 + samtools + bcftools
 - Re-aligns consensus to the full DB and reports whether the consensus
   best-hit matches the original reference and the percent identity.

Notes:
 - Requires minimap2, samtools and bcftools in PATH (or use conda run -p ./env ...)
 - Assumes PAF is in minimap2 format (nmatch in column 10, mapq in column 12 index 11)
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
    with opener(reads_gz, 'rt') as inf:
        first = None
        while True:
            pos = inf.tell()
            line = inf.readline()
            if not line:
                break
            if line.strip() == '':
                continue
            first = line[0]
            inf.seek(pos)
            break

        inf.close()

    if first == '>':
        with opener(reads_gz, 'rt') as inf, open(out_fastq, 'wt') as outf:
            cur_name = None
            cur_seq_parts = []
            for line in inf:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith('>'):
                    if cur_name is not None:
                        if cur_name in names:
                            outf.write(f">{cur_name}\n")
                            outf.write('\n'.join(cur_seq_parts) + '\n')
                            found += 1
                            names.remove(cur_name)
                            if not names:
                                break
                    cur_name = line[1:].split()[0]
                    cur_seq_parts = []
                else:
                    cur_seq_parts.append(line)
            if cur_name is not None and cur_name in names:
                outf.write(f">{cur_name}\n")
                outf.write('\n'.join(cur_seq_parts) + '\n')
                found += 1
    else:
        opener = gzip.open if str(reads_gz).endswith('.gz') else open
        with opener(reads_gz, 'rt') as inf, open(out_fastq, 'wt') as outf:
            while True:
                hdr = inf.readline()
                if not hdr:
                    break
                seq = inf.readline()
                plus = inf.readline()
                qual = inf.readline()
                if not qual:
                    break
                token = hdr.strip().split()[0]
                if token.startswith('@'):
                    token = token[1:]
                if token in names:
                    outf.write(hdr)
                    outf.write(seq)
                    outf.write(plus)
                    outf.write(qual)
                    found += 1
                    names.remove(token)
                if not names:
                    break

    return found


def safe_name(name: str) -> str:
    import re
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def make_single_ref_fasta(ref_id, seq, out_path):
    with open(out_path, "wt") as fh:
        fh.write(f">{ref_id}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")


def build_consensus(ref_fa, reads_fastq, work_dir, threads=4):
    ref_fa = str(ref_fa)
    reads_fastq = str(reads_fastq)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    sam = work_dir / "aln.sam"
    bam = work_dir / "aln.bam"
    sorted_bam = work_dir / "aln.sorted.bam"
    vcf = work_dir / "calls.vcf.gz"
    consensus = work_dir / "consensus.fa"

    cmd = ["minimap2", "-a", "-x", "map-ont", "-t", str(threads), ref_fa, reads_fastq]
    with open(sam, "wb") as samfh:
        print("Running:", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, stdout=samfh)
        except subprocess.CalledProcessError:
            print("ERROR: minimap2 failed", file=sys.stderr)
            return None

    run_cmd(["samtools", "view", "-bS", str(sam), "-o", str(bam)])
    run_cmd(["samtools", "sort", "-o", str(sorted_bam), str(bam)])
    run_cmd(["samtools", "index", str(sorted_bam)])

    mpileup_cmd = f"bcftools mpileup -f {ref_fa} {str(sorted_bam)} | bcftools call -mv -Oz -o {str(vcf)} -"
    run_cmd(mpileup_cmd, shell=True)
    run_cmd(["bcftools", "index", str(vcf)])

    run_cmd(["bcftools", "consensus", "-f", ref_fa, str(vcf), "-o", str(consensus)])

    return consensus


def parse_paf_best_hit(paf_file):
    best = None
    with open(paf_file, "rt") as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            qname = cols[0]
            tname = cols[5]
            try:
                nmatch = int(cols[9])
            except ValueError:
                nmatch = 0
            try:
                aln_len = int(cols[10])
            except ValueError:
                aln_len = 0
            if best is None or nmatch > best[1]:
                best = (tname, nmatch, aln_len)
    return best


def main():
    parser = argparse.ArgumentParser(description="Consensus-based validation for each reference")
    parser.add_argument("--paf_file", required=True)
    parser.add_argument("--reads_file", required=True)
    parser.add_argument("--db_file", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--mapq", type=int, default=40, help="Minimum MAPQ to consider an alignment")
    parser.add_argument("--min_reads", type=int, default=20, help="Minimum reads per target to attempt consensus")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    paf = Path(args.paf_file)
    reads = Path(args.reads_file)
    db = Path(args.db_file)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Loading reference DB...")
    refs = load_fasta(str(db))
    print(f"Loaded {len(refs)} reference sequences")

    print("Grouping reads by reference from PAF (MAPQ>={})...".format(args.mapq))
    groups = parse_paf_group_reads(str(paf), mapq_threshold=args.mapq)

    report_rows = []
    for ref_id, read_names in groups.items():
        read_count = len(read_names)
        print(f"Reference {ref_id}: {read_count} high-MAPQ reads")
        if read_count < args.min_reads:
            print(f"  Skipping {ref_id}: not enough reads (need {args.min_reads})")
            report_rows.append((ref_id, read_count, 0, "", 0.0, "TOO_FEW_READS"))
            continue

        safe = safe_name(ref_id)
        reads_fastq = outdir / f"{safe}.reads.fastq"
        print(f"  Extracting reads to {reads_fastq} ...")
        found = extract_reads_by_names(str(reads), read_names, str(reads_fastq))
        print(f"  Extracted {found}/{read_count} reads")

        if ref_id not in refs:
            print(f"  WARNING: reference {ref_id} not found in DB; skipping", file=sys.stderr)
            report_rows.append((ref_id, found, 0, "", 0.0, "NO_REF"))
            continue
        single_ref_fa = outdir / f"{safe}.ref.fa"
        make_single_ref_fasta(ref_id, refs[ref_id], single_ref_fa)

        work = outdir / f"work_{safe}"
        print(f"  Building consensus for {ref_id} ...")
        consensus_fa = build_consensus(single_ref_fa, reads_fastq, work, threads=args.threads)
        if consensus_fa is None or not Path(consensus_fa).exists():
            print(f"  Consensus build failed for {ref_id}")
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
        except subprocess.CalledProcessError:
            print(f"  Error: minimap2 failed aligning consensus for {ref_id}", file=sys.stderr)
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
