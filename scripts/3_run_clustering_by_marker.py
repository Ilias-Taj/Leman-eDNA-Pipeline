#!/usr/bin/env python3
"""
3_run_clustering_by_marker.py - CLUSTERING + CHIMERA DETECTION

Performs VSEARCH clustering separately for 18S and COI markers at 95% identity.
Generates consensus sequences for error correction and removes PCR chimeras.

Prerequisite: Run 2_classify_markers.py first

Usage:
    python scripts/3_run_clustering_by_marker.py \
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \
        --output_dir out/Water_eDNA_18S_COI_14_01_26 \
        --identity 0.95 \
        --threads 12

Output:
    Creates temp_clustering/ folder with:
    - consensus_18S_clean.fasta (chimera-filtered consensus sequences)
    - consensus_COI_clean.fasta (chimera-filtered consensus sequences)
    - global_otu_assignment.txt (read-to-OTU mapping for 18S)
    - global_otu_assignment_COI.txt (read-to-OTU mapping for COI)
"""

import argparse
import sys
import subprocess
import shutil
from pathlib import Path
import gzip
import os

def remove_chimeras(working_dir, marker, vsearch_path):
    """
    Detect and remove chimeric OTUs using VSEARCH --uchime_denovo.
    """
    consensus_in = working_dir / f"consensus_{marker}.fasta"
    # Output files
    consensus_clean = working_dir / f"consensus_{marker}_clean.fasta"
    chimeras_file = working_dir / f"chimeras_{marker}.fasta"
    
    cmd = [
        vsearch_path,
        "--uchime_denovo", str(consensus_in),
        "--nonchimeras", str(consensus_clean),
        "--chimeras", str(chimeras_file),
        "--xsize"  # Removes abundance annotations for processing
    ]
    
    print(f"Running Chimera Detection for {marker}...")
    subprocess.run(cmd, check=True)
    
    # IMPORTANT: You must now update your OTU table to remove these IDs
    # or simply use the 'consensus_clean.fasta' for taxonomy assignment later.
    return consensus_clean, chimeras_file

def load_chimera_ids(chimeras_fasta):
    """Load chimera centroid IDs from a FASTA file."""
    chimera_ids = set()
    if not chimeras_fasta.exists():
        return chimera_ids
    with open(chimeras_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip().split()[0]
                clean_id = header.split(';')[0]
                chimera_ids.add(clean_id)
    return chimera_ids

def check_vsearch():
    """Check if vsearch is installed."""
    # Try multiple locations
    candidates = [
        Path("./env/bin/vsearch"),
        Path("/opt/homebrew/bin/vsearch"),
        Path("/usr/local/bin/vsearch"),
    ]
    
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    
    vsearch = shutil.which("vsearch")
    if vsearch:
        return vsearch
    
    print("ERROR: 'vsearch' tool not found. Please install it:", file=sys.stderr)
    print("  brew install vsearch  (or conda install -c bioconda vsearch)", file=sys.stderr)
    sys.exit(1)

def concatenate_reads_by_marker(input_dir, marker, output_fasta):
    """
    Find all filtered_reads_{marker}.fastq.gz, convert to FASTA, 
    and merge into one file with barcode tags in headers.
    """
    print(f"Collecting {marker} reads from {input_dir}...")
    
    total_reads = 0
    with open(output_fasta, 'w') as out_f:
        # Sort to ensure consistent order
        for sample_dir in sorted(input_dir.iterdir()):
            if not sample_dir.is_dir() or sample_dir.name in ("logs", "validation", "merged", "temp_clustering"):
                continue
            
            fq_path = sample_dir / f"filtered_reads_{marker}.fastq.gz"
            if not fq_path.exists():
                continue
            
            barcode = sample_dir.name
            
            # Open gzip file and process FASTQ (4 lines per read)
            with gzip.open(fq_path, 'rt') as in_f:
                lines = []
                for line in in_f:
                    lines.append(line.rstrip('\n'))
                    
                    # FASTQ format: @header, sequence, +, quality
                    if len(lines) == 4:
                        header = lines[0][1:].split()[0]  # Remove @ and get first part
                        seq = lines[1]
                        
                        # Write as FASTA: >readname|barcode|marker
                        out_f.write(f">{header}|{barcode}|{marker}\n{seq}\n")
                        total_reads += 1
                        lines = []
    
    print(f"Concatenated {total_reads} {marker} reads into {output_fasta}")
    return total_reads

def run_vsearch_clustering(input_fasta, working_dir, marker, identity, threads, vsearch_path):
    """Run VSEARCH --cluster_fast for a specific marker."""
    consensus_file = working_dir / f"consensus_{marker}.fasta"
    uc_file = working_dir / f"clusters_{marker}.uc"
    
    cmd = [
        vsearch_path,
        "--cluster_fast", str(input_fasta),
        "--id", str(identity),
        "--minseqlength", "300",
        "--maxseqlength", "3000",
        "--consout", str(consensus_file),
        "--uc", str(uc_file),
        "--threads", str(threads),
        "--qmask", "none",
        "--sizeout"
    ]
    
    print(f"Running VSEARCH for {marker}: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"VSEARCH Failed for {marker}:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)
        
    print(f"VSEARCH clustering complete for {marker}.")
    return uc_file

def parse_uc_to_assignment(uc_file, output_file, marker, chimera_ids=None):
    """
    Parse UC output to global OTU assignments.
    Format: read_name<TAB>otu_id<TAB>barcode<TAB>marker
    
    Read ID format from concatenation: {uuid}|{barcode}|{marker}
    """
    otu_mapping = {}
    otu_counter = 1
    written_count = 0
    chimera_ids = chimera_ids or set()
    
    with open(uc_file, 'r') as in_f:
        for line in in_f:
            cols = line.strip().split('\t')
            
            if not cols:
                continue
            
            record_type = cols[0]
            
            if record_type == 'S':
                # Cluster seed (centroid)
                read_full_id = cols[8]
                if read_full_id in chimera_ids:
                    continue
                otu_id = f"OTU_{marker}_{otu_counter:06d}"
                otu_mapping[read_full_id] = otu_id
                otu_counter += 1
            
            elif record_type == 'H':
                # Hit (member of cluster)
                read_full_id = cols[8]
                centroid_full_id = cols[9]
                if centroid_full_id in chimera_ids:
                    continue
                otu_id = otu_mapping.get(centroid_full_id)
                
                if otu_id and read_full_id not in otu_mapping:
                    otu_mapping[read_full_id] = otu_id
    
    # Write assignments file with proper parsing
    with open(output_file, 'w') as f_out:
        for read_full_id, clean_otu_id in otu_mapping.items():
            # Parse: {uuid}|{barcode}|{marker}
            parts = read_full_id.split('|')
            if len(parts) >= 3:
                uuid = parts[0]
                barcode = parts[1]
                m = parts[2]
                f_out.write(f"{uuid}\t{clean_otu_id}\t{barcode}\t{m}\n")
                written_count += 1
    
    print(f"Generated assignments for {otu_counter-1} {marker} OTUs, wrote {written_count} assignments.")
    return otu_counter - 1

def main():
    parser = argparse.ArgumentParser(description="Global OTU clustering by marker using VSEARCH")
    parser.add_argument("--input_dir", required=True, help="Input directory containing barcodes")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--identity", type=float, default=0.95, help="Identity threshold")
    parser.add_argument("--threads", type=int, default=4, help="Threads for vsearch")
    
    args = parser.parse_args()
    
    vsearch_path = check_vsearch()
    
    input_dir = Path(args.input_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Temp files
    temp_dir = out_dir / "temp_clustering"
    temp_dir.mkdir(exist_ok=True)
    
    # Process each marker separately
    all_otus = {}
    
    for marker in ["18S", "COI"]:
        print(f"\n{'='*60}")
        print(f"PROCESSING MARKER: {marker}")
        print(f"{'='*60}")
        
        merged_fasta = temp_dir / f"all_reads_{marker}.fasta"
        
        # 1. Merge Reads
        num_reads = concatenate_reads_by_marker(input_dir, marker, merged_fasta)
        if num_reads == 0:
            print(f"[SKIP] No {marker} reads found")
            continue
        
        # 2. Run VSEARCH
        uc_file = run_vsearch_clustering(merged_fasta, temp_dir, marker, args.identity, args.threads, vsearch_path)

        # 3. Chimera detection on centroids
        _, chimeras_file = remove_chimeras(temp_dir, marker, vsearch_path)
        chimera_ids = load_chimera_ids(chimeras_file)
        
        # 4. Parse Output (skip chimeric OTUs)
        final_assignment = out_dir / f"global_otu_assignment_{marker}.txt"
        num_otus = parse_uc_to_assignment(uc_file, final_assignment, marker, chimera_ids=chimera_ids)
        all_otus[marker] = num_otus
        
        print(f"✓ {marker} OTU assignment file: {final_assignment}")
    
    # Merge all marker assignments into one file
    print(f"\n{'='*60}")
    print("MERGING MARKER ASSIGNMENTS")
    print(f"{'='*60}")
    
    merged_assignment = out_dir / "global_otu_assignment.txt"
    with open(merged_assignment, 'w') as out_f:
        # Header
        out_f.write("read_name\totu_id\tbarcode\tmarker\n")
        
        # Append each marker's assignments
        for marker in ["18S", "COI"]:
            assignment_file = out_dir / f"global_otu_assignment_{marker}.txt"
            if assignment_file.exists():
                with open(assignment_file, 'r') as in_f:
                    for line in in_f:
                        out_f.write(line)
    
    # Print summary
    print(f"\nPipeline Finished!")
    print(f"18S OTUs: {all_otus.get('18S', 0)}")
    print(f"COI OTUs: {all_otus.get('COI', 0)}")
    print(f"Total OTUs: {sum(all_otus.values())}")
    print(f"Global assignment file (with markers): {merged_assignment}")

if __name__ == "__main__":
    main()
