#!/usr/bin/env python3
"""
6_blast_top_otus.py - BLAST TOP OTUs AGAINST NCBI (ROBUST VERSION)

Extracts the top N most abundant OTUs from the abundance matrix
and runs remote BLAST searches against NCBI GenBank for taxonomy validation.

Robustly handles ID mapping from OTU names to consensus sequences by
checking multiple candidate reads until finding one that matches the FASTA file.

Prerequisite: Run 4_merge_otu_tables_by_marker.py first

Usage:
    # BLAST top 10 COI OTUs
    python scripts/6_blast_top_otus.py \
        --matrix out/Water_eDNA_18S_COI_14_01_26/merged/otu_relative_abundance_COI.csv \
        --fasta out/Water_eDNA_18S_COI_14_01_26/temp_clustering/consensus_COI_clean.fasta \
        --otu_assignment out/Water_eDNA_18S_COI_14_01_26/global_otu_assignment_COI.txt \
        --marker COI \
        --top_n 10
    
    # BLAST top 10 18S OTUs
    python scripts/6_blast_top_otus.py \
        --matrix out/Water_eDNA_18S_COI_14_01_26/merged/otu_relative_abundance_18S.csv \
        --fasta out/Water_eDNA_18S_COI_14_01_26/temp_clustering/consensus_18S_clean.fasta \
        --otu_assignment out/Water_eDNA_18S_COI_14_01_26/global_otu_assignment.txt \
        --marker 18S \
        --top_n 10

Output:
    Results are saved to {output_dir}/blast_top{N}_{marker}.txt
    Default output_dir: {input_dir}/blast_results/

Note:
    BLAST queries are rate-limited (3 seconds between queries) to comply
    with NCBI usage policies. Expect ~1-2 minutes per sequence.
"""

import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import argparse
import sys
import time
from pathlib import Path

def load_otu_mapping_robust(otu_assignment_file):
    """
    Maps OTU_ID -> List of candidate Read IDs.
    
    Stores multiple read IDs per OTU to find the centroid that exists
    in the actual FASTA file (which only contains consensus sequences).
    """
    otu_to_reads = {}
    print("  Loading assignment file (this may take a moment)...")
    
    with open(otu_assignment_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                read_full_id = parts[0]  # uuid|barcode|marker
                otu_id = parts[1]
                
                if otu_id not in otu_to_reads:
                    otu_to_reads[otu_id] = set()
                
                # Store first 5 candidates per OTU (optimization to save RAM)
                if len(otu_to_reads[otu_id]) < 5:
                    otu_to_reads[otu_id].add(read_full_id)
    
    return otu_to_reads

def get_sequence_for_otu(fasta_path, otu_id, otu_to_reads_map):
    """
    Scans FASTA file to find a sequence matching one of the reads assigned to this OTU.
    Returns (sequence, found_id) or (None, None) if not found.
    """
    if otu_id not in otu_to_reads_map:
        return None, None
    
    candidates = otu_to_reads_map[otu_id]
    
    with open(fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # VSEARCH consensus headers look like: >centroid=uuid|barcode|marker;size=100
            header_clean = record.id.split(';')[0].replace('centroid=', '')
            
            # Try matching different parts:
            # 1. Full ID match (uuid|barcode|marker)
            if header_clean in candidates:
                return str(record.seq), header_clean
            
            # 2. Just the UUID part (before first |)
            uuid_part = header_clean.split('|')[0]
            if uuid_part in candidates:
                return str(record.seq), uuid_part
            
            # 3. Try any candidate that starts with this UUID
            for candidate in candidates:
                if header_clean.startswith(candidate) or candidate.startswith(uuid_part):
                    return str(record.seq), header_clean
    
    return None, None

def main():
    parser = argparse.ArgumentParser(
        description="BLAST top N OTUs against NCBI GenBank for taxonomy validation",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--matrix", required=True, 
                        help="Path to relative abundance matrix CSV")
    parser.add_argument("--fasta", required=True, 
                        help="Path to consensus FASTA file (consensus_*_clean.fasta)")
    parser.add_argument("--otu_assignment", required=True, 
                        help="Path to global OTU assignment file for OTU-to-centroid mapping")
    parser.add_argument("--marker", required=True, 
                        help="Marker name (18S or COI)")
    parser.add_argument("--top_n", type=int, default=10, 
                        help="Number of top abundant OTUs to BLAST (default: 10)")
    parser.add_argument("--output_dir", default=None, 
                        help="Output directory for results (default: {input_dir}/blast_results/)")
    args = parser.parse_args()

    # Set up output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = Path(args.matrix).parent.parent / "blast_results"
    output_dir.mkdir(exist_ok=True, parents=True)
    
    output_file = output_dir / f"blast_top{args.top_n}_{args.marker}.txt"

    print(f"--- BLASTing Top {args.top_n} {args.marker} OTUs ---")
    print(f"Results will be saved to: {output_file}")

    # 1. Load robust OTU to reads mapping
    print("Loading OTU to reads mapping...")
    otu_to_reads = load_otu_mapping_robust(args.otu_assignment)
    print(f"  Loaded {len(otu_to_reads)} OTU mappings")

    # 2. Load Matrix and find Top OTUs
    try:
        df = pd.read_csv(args.matrix, index_col=0)
        # Sum all sample columns
        df['total_reads'] = df.sum(axis=1)
        
        # Sort and take top N
        top_otus = df.sort_values('total_reads', ascending=False).head(args.top_n)
        top_ids = list(top_otus.index)
        
        print(f"Identified Top {args.top_n} OTUs (from {len(df)} total).")
    except Exception as e:
        print(f"Error loading matrix: {e}", file=sys.stderr)
        sys.exit(1)

    # 3. Extract sequences from FASTA for top OTUs (only ~10, no indexing needed)

    sequences = {}
    for otu_id in top_ids:
        seq, found_id = get_sequence_for_otu(args.fasta, otu_id, otu_to_reads)
        if seq:
            sequences[otu_id] = seq
            print(f"  > Loaded sequence for {otu_id} ({len(seq)} bp from {found_id[:20]}...)")
        else:
            print(f"  [WARN] Could not find sequence for {otu_id}", file=sys.stderr)

    # 4. Run Remote BLAST
    print("\nStarting Remote BLAST (This takes 1-2 minutes per sequence)...")
    
    # Open output file for writing
    with open(output_file, 'w') as out_f:
        # Write header
        header = f"{'OTU ID':<20} | {'Reads':<10} | {'Top Hit Species':<40} | {'Identity':<10} | {'E-Value'}\n"
        separator = "-" * 110 + "\n"
        
        print(header + separator, end="")
        out_f.write(f"=== BLAST Results for Top {args.top_n} {args.marker} OTUs ===\n\n")
        out_f.write(header + separator)

        for otu_id in top_ids:
            if otu_id not in top_otus.index:
                continue
                
            count = top_otus.loc[otu_id, 'total_reads']
            
            if otu_id not in sequences:
                result_line = f"{otu_id:<20} | {count:<10.4f} | SEQUENCE NOT FOUND                      | -          | -\n"
                print(result_line, end="")
                out_f.write(result_line)
                continue
                
            seq = sequences[otu_id]
            
            try:
                # BLASTN against 'nt' (Nucleotide collection)
                print(f"  BLASTing {otu_id}...", end=" ", flush=True)
                result_handle = NCBIWWW.qblast("blastn", "nt", seq, hitlist_size=1)
                blast_record = NCBIXML.read(result_handle)
                
                if blast_record.alignments:
                    alignment = blast_record.alignments[0]
                    hsp = alignment.hsps[0]
                    species_name = alignment.title.split("|")[-1].strip()
                    # Clean up species name (often contains extensive info)
                    if " " in species_name:
                        # heuristic to grab just the species name part
                        parts = species_name.split()
                        if len(parts) > 2:
                            species_name = " ".join(parts[0:3])
                    
                    ident = f"{100 * hsp.identities / hsp.align_length:.1f}%"
                    print(f"done")
                    result_line = f"{otu_id:<20} | {count:<10.4f} | {species_name:<40} | {ident:<10} | {hsp.expect:.2e}\n"
                    print(result_line, end="")
                    out_f.write(result_line)
                else:
                    print(f"done")
                    result_line = f"{otu_id:<20} | {count:<10.4f} | NO MATCH FOUND                          | -          | -\n"
                    print(result_line, end="")
                    out_f.write(result_line)
                
                # Be polite to NCBI servers
                time.sleep(3)
                
            except Exception as e:
                print(f"error")
                result_line = f"{otu_id:<20} | ERROR: {e}\n"
                print(result_line, end="")
                out_f.write(result_line)
    
    print(f"\n[OK] Results saved to: {output_file}")

if __name__ == "__main__":
    main()
