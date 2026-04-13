#!/usr/bin/env python3
"""
7_comprehensive_taxonomy_summary.py - COMPREHENSIVE TAXONOMY REPORT

Combines local taxonomy (SILVA/SINTAX), optional NCBI BLAST results,
and abundance data into a single CSV file for downstream analysis.

Prerequisite: Run 5_assign_taxonomy.py first

Usage:
    # Generate comprehensive CSV with BLAST validation for top 50 OTUs
    python scripts/7_comprehensive_taxonomy_summary.py \
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \
        --blast_n 50
    
    # Generate CSV from existing SINTAX taxonomy only (no BLAST)
    python scripts/7_comprehensive_taxonomy_summary.py \
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \
        --skip_blast

Output:
    Saves comprehensive taxonomy CSV files to:
    - {input_dir}/comprehensive_taxonomy_18S.csv
    - {input_dir}/comprehensive_taxonomy_COI.csv
"""

import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import argparse
import sys
import time
from pathlib import Path
import re

def load_otu_to_centroid_mapping(otu_assignment_file):
    """Load mapping from OTU IDs to centroid IDs."""
    otu_to_centroid = {}
    
    with open(otu_assignment_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('read_name'):
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    centroid_id = parts[0]
                    otu_id = parts[1]
                    if otu_id not in otu_to_centroid:
                        otu_to_centroid[otu_id] = centroid_id
    
    return otu_to_centroid

def clean_taxon_name(taxon):
    """
    Clean MIDORI2/SILVA taxon names:
    1. Strip rank prefixes (e.g., "order_Vannellidae_95227" -> "Vannellidae_95227")
    2. Strip trailing NCBI taxon IDs (e.g., "Arthropoda_6656" -> "Arthropoda")
    3. Strip trailing underscores (e.g., "Metazoa_" -> "Metazoa")
    """
    if not taxon:
        return taxon
    # Strip rank prefixes
    taxon = re.sub(r'^(kingdom|phylum|class|order|family|genus|species)_', '', taxon, flags=re.IGNORECASE)
    # Strip trailing NCBI taxon IDs (digits after last underscore)
    parts = taxon.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        taxon = parts[0]
    # Strip trailing underscores
    taxon = taxon.rstrip('_')
    return taxon

def parse_silva_taxonomy(taxonomy_file):
    """Parse SILVA/MIDORI2 taxonomy assignments from SINTAX output."""
    taxonomy_dict = {}
    
    with open(taxonomy_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            
            # Extract centroid ID from header
            # Handles both formats:
            #   centroid=UUID|barcode|marker  (old VSEARCH --consout)
            #   UUID|barcode|marker           (current pipeline)
            header = parts[0].split(';')[0]  # strip ;size= if present
            match = re.search(r'centroid=([a-f0-9\-]+)', header)
            if match:
                centroid_id = match.group(1)
            else:
                # Use first field (UUID) from pipe-delimited header
                centroid_id = header.split('|')[0].strip()
            full_taxonomy = parts[1] if len(parts) > 1 else ""
            
            # Parse taxonomy levels
            tax_dict = {'domain': '', 'phylum': '', 'class': '', 'order': '', 
                       'family': '', 'genus': '', 'species': ''}
            
            for match in re.finditer(r'([dkpcofgs]):([^,()]+)', full_taxonomy):
                level_code = match.group(1)
                taxon = match.group(2).strip()
                taxon = clean_taxon_name(taxon)
                
                level_map = {'d': 'domain', 'k': 'kingdom', 'p': 'phylum', 
                           'c': 'class', 'o': 'order', 'f': 'family', 
                           'g': 'genus', 's': 'species'}
                
                if level_code in level_map:
                    tax_dict[level_map[level_code]] = taxon
            
            # MIDORI2 (COI) uses k: (kingdom) instead of d: (domain).
            # Fall back to kingdom when domain is empty.
            if not tax_dict.get('domain') and tax_dict.get('kingdom'):
                tax_dict['domain'] = tax_dict['kingdom']
            
            taxonomy_dict[centroid_id] = tax_dict
    
    return taxonomy_dict

def run_blast_batch(sequences, otu_ids):
    """Run BLAST for a batch of sequences."""
    blast_results = {}
    
    print(f"\nBLASTing {len(sequences)} sequences...")
    
    for i, otu_id in enumerate(otu_ids, 1):
        if otu_id not in sequences:
            blast_results[otu_id] = {'species': 'N/A', 'identity': 0, 'evalue': 'N/A'}
            continue
        
        seq = sequences[otu_id]
        
        try:
            print(f"  [{i}/{len(otu_ids)}] BLASTing {otu_id}...", end=" ", flush=True)
            result_handle = NCBIWWW.qblast("blastn", "nt", seq, hitlist_size=1)
            blast_record = NCBIXML.read(result_handle)
            
            if blast_record.alignments:
                alignment = blast_record.alignments[0]
                hsp = alignment.hsps[0]
                species_name = alignment.title.split("|")[-1].strip()
                
                # Clean species name
                if " " in species_name:
                    parts = species_name.split()
                    if len(parts) > 2:
                        species_name = " ".join(parts[0:3])
                
                identity = 100 * hsp.identities / hsp.align_length
                
                blast_results[otu_id] = {
                    'species': species_name,
                    'identity': round(identity, 1),
                    'evalue': f"{hsp.expect:.2e}"
                }
                print(f"done ({identity:.1f}%)")
            else:
                blast_results[otu_id] = {'species': 'No match', 'identity': 0, 'evalue': 'N/A'}
                print("no match")
            
            time.sleep(3)  # Be polite to NCBI
            
        except Exception as e:
            print(f"error: {e}")
            blast_results[otu_id] = {'species': f'Error: {str(e)[:50]}', 'identity': 0, 'evalue': 'N/A'}
    
    return blast_results

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="Input directory")
    parser.add_argument("--blast_n", type=int, default=0, help="Number of top OTUs to BLAST per marker (0=skip)")
    parser.add_argument("--skip_blast", action='store_true', help="Skip BLAST, use only local taxonomy")
    parser.add_argument("--markers", default=None,
                        help="Comma-separated list of markers (default: auto-detect from merged/ files)")
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    output_dir = input_dir / "taxonomy_summary"
    output_dir.mkdir(exist_ok=True)
    
    # Determine markers to process
    if args.markers:
        markers_to_process = [m.strip().upper() for m in args.markers.split(",")]
    else:
        # Auto-detect from existing abundance files
        merged_dir = input_dir / "merged"
        markers_to_process = []
        for candidate in ["18S", "COI", "JEDI"]:
            if (merged_dir / f"otu_relative_abundance_{candidate}.csv").exists():
                markers_to_process.append(candidate)
        if not markers_to_process:
            markers_to_process = ["18S", "COI"]
    
    print("=" * 80)
    print("COMPREHENSIVE TAXONOMY SUMMARY")
    print(f"Markers: {', '.join(markers_to_process)}")
    print("=" * 80)
    
    for marker in markers_to_process:
        print(f"\n{'=' * 80}")
        print(f"PROCESSING {marker}")
        print(f"{'=' * 80}")
        
        # File paths
        abundance_file = input_dir / f"merged/otu_relative_abundance_{marker}.csv"
        taxonomy_file = input_dir / f"taxonomy/taxonomy_{marker}.txt"
        consensus_file = input_dir / f"temp_clustering/consensus_{marker}_clean.fasta"
        otu_assignment_file = input_dir / f"global_otu_assignment_{marker}.txt"
        
        # Check files exist
        if not abundance_file.exists():
            print(f"  ⚠ Abundance file not found: {abundance_file}")
            continue
        
        # 1. Load abundance data
        print("\n[1/5] Loading abundance data...")
        abundance_df = pd.read_csv(abundance_file, index_col=0)
        abundance_df['total_abundance'] = abundance_df.sum(axis=1)
        abundance_df = abundance_df.sort_values('total_abundance', ascending=False)
        print(f"  Loaded {len(abundance_df)} OTUs")
        
        # 2. Load OTU to centroid mapping
        print("[2/5] Loading OTU to centroid mapping...")
        otu_to_centroid = load_otu_to_centroid_mapping(otu_assignment_file)
        print(f"  Loaded {len(otu_to_centroid)} mappings")
        
        # 3. Load SILVA taxonomy
        print("[3/5] Loading local taxonomy assignments...")
        silva_taxonomy = {}
        if taxonomy_file.exists():
            silva_taxonomy = parse_silva_taxonomy(taxonomy_file)
            print(f"  Loaded {len(silva_taxonomy)} taxonomy assignments")
        else:
            print(f"  ⚠ Taxonomy file not found")
        
        # 4. Run BLAST if requested
        blast_results = {}
        if not args.skip_blast and args.blast_n > 0:
            print(f"[4/5] Running BLAST on top {args.blast_n} OTUs...")
            
            # Get top N OTUs
            top_otus = abundance_df.head(args.blast_n).index.tolist()
            
            # Load sequences
            sequences = {}
            with open(consensus_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    # Strip ;size= and centroid= prefixes, then take UUID
                    raw_id = record.id.split(';')[0]
                    centroid_id = raw_id.split('|')[0].replace('centroid=', '')
                    
                    for otu_id in top_otus:
                        if otu_id in otu_to_centroid and otu_to_centroid[otu_id] == centroid_id:
                            sequences[otu_id] = str(record.seq)
                            break
            
            print(f"  Loaded {len(sequences)} sequences")
            blast_results = run_blast_batch(sequences, top_otus)
        else:
            print("[4/5] Skipping BLAST")
        
        # 5. Create comprehensive summary
        print("[5/5] Creating summary CSV...")
        
        summary_rows = []
        for otu_id in abundance_df.index:
            row = {
                'OTU_ID': otu_id,
                'Total_Abundance': abundance_df.loc[otu_id, 'total_abundance'],
                'Rank': len(summary_rows) + 1
            }
            
            # Add per-sample abundances
            for sample in abundance_df.columns:
                if sample != 'total_abundance':
                    row[f'Sample_{sample}'] = abundance_df.loc[otu_id, sample]
            
            # Add SILVA taxonomy
            centroid_id = otu_to_centroid.get(otu_id, '')
            if centroid_id in silva_taxonomy:
                tax = silva_taxonomy[centroid_id]
                row['SILVA_Domain'] = tax.get('domain', '')
                row['SILVA_Phylum'] = tax.get('phylum', '')
                row['SILVA_Class'] = tax.get('class', '')
                row['SILVA_Order'] = tax.get('order', '')
                row['SILVA_Family'] = tax.get('family', '')
                row['SILVA_Genus'] = tax.get('genus', '')
                row['SILVA_Species'] = tax.get('species', '')
            else:
                row['SILVA_Domain'] = ''
                row['SILVA_Phylum'] = ''
                row['SILVA_Class'] = ''
                row['SILVA_Order'] = ''
                row['SILVA_Family'] = ''
                row['SILVA_Genus'] = ''
                row['SILVA_Species'] = ''
            
            # Add BLAST results if available
            if otu_id in blast_results:
                row['NCBI_TopHit'] = blast_results[otu_id]['species']
                row['NCBI_Identity'] = blast_results[otu_id]['identity']
                row['NCBI_Evalue'] = blast_results[otu_id]['evalue']
            else:
                row['NCBI_TopHit'] = ''
                row['NCBI_Identity'] = ''
                row['NCBI_Evalue'] = ''
            
            summary_rows.append(row)
        
        summary_df = pd.DataFrame(summary_rows)
        
        # Save to CSV
        output_file = output_dir / f"comprehensive_taxonomy_{marker}.csv"
        summary_df.to_csv(output_file, index=False)
        print(f"\n✓ Saved: {output_file}")
        print(f"  {len(summary_df)} OTUs with taxonomy and abundance data")
        
        # Print summary statistics
        print(f"\nSummary Statistics for {marker}:")
        print(f"  Total OTUs: {len(summary_df)}")
        assigned = summary_df[summary_df['SILVA_Phylum'] != ''].shape[0]
        print(f"  With SILVA taxonomy: {assigned} ({100*assigned/len(summary_df):.1f}%)")
        if blast_results:
            blasted = summary_df[summary_df['NCBI_TopHit'] != ''].shape[0]
            print(f"  With NCBI BLAST: {blasted} ({100*blasted/len(summary_df):.1f}%)")
    
    print("\n" + "=" * 80)
    print("COMPREHENSIVE TAXONOMY SUMMARY COMPLETE")
    print(f"Output directory: {output_dir}")
    print("=" * 80)

if __name__ == "__main__":
    main()
