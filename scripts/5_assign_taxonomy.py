#!/usr/bin/env python3
"""
5_assign_taxonomy.py - ASSIGN TAXONOMY TO CONSENSUS OTU SEQUENCES

Uses VSEARCH --sintax to assign taxonomy to clean consensus sequences.
Requires reference databases in VSEARCH .udb format.

Prerequisites:
    1. Run 3_run_clustering_by_marker.py first
    2. Download and prepare reference databases:
       - 18S: SILVA database (https://www.drive5.com/usearch/manual/sintax_downloads.html)
       - COI: Database to be determined 
    3. Convert to .udb format:
       wget https://www.drive5.com/sintax/silva_18s_v123.fa.gz
       gunzip silva_18s_v123.fa.gz
       vsearch --makeudb_usearch silva_18s_v123.fa --output silva_18s.udb

Usage Examples:
    # With both databases
    python scripts/5_assign_taxonomy.py \
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \
        --db_18S refs/silva_18s_v123.udb \
        --db_COI refs/midori2_COI.udb \
        --threads 12
    
    # 18S only
    python scripts/5_assign_taxonomy.py \
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \
        --db_18S refs/silva_18s_v123.udb \
        --threads 12
    
    # COI only
    python scripts/5_assign_taxonomy.py \
        --input_dir out/Water_eDNA_18S_COI_14_01_26 \
        --db_COI refs/midori2_COI.udb \
        --threads 12

Output:
    Saves taxonomy assignments to:
    - {input_dir}/taxonomy/taxonomy_18S.txt
    - {input_dir}/taxonomy/taxonomy_COI.txt
"""

import argparse
import sys
import subprocess
from pathlib import Path

# Ensure scripts/ is on the import path so utils.py can be found from any working directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from db_tag import label_from_path
from utils import find_tool

def run_sintax(consensus_fasta, output_file, db_path, vsearch_path, threads, confidence=0.8):
    """
    Run VSEARCH --sintax for taxonomy assignment.
    
    Args:
        consensus_fasta: Path to clean consensus sequences
        output_file: Output taxonomy file
        db_path: Path to reference database (.udb file)
        vsearch_path: Path to vsearch executable
        threads: Number of threads
        confidence: Confidence threshold (default 0.8)
    """
    if not Path(db_path).exists():
        print(f"WARNING: Database not found: {db_path}", file=sys.stderr)
        print(f"Skipping taxonomy assignment. Please download appropriate database.", file=sys.stderr)
        return False
    
    cmd = [
        vsearch_path,
        "--sintax", str(consensus_fasta),
        "--db", str(db_path),
        "--tabbedout", str(output_file),
        "--sintax_cutoff", str(confidence),
        "--threads", str(threads)
    ]
    
    print(f"Running SINTAX taxonomy assignment...")
    print(f"  Input: {consensus_fasta}")
    print(f"  Database: {db_path}")
    print(f"  Output: {output_file}")
    print(f"  Confidence threshold: {confidence}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"SINTAX Failed:\n{result.stderr}", file=sys.stderr)
        return False
    
    print(f"[OK] Taxonomy assignment complete.")
    return True

def parse_sintax_output(sintax_file, marker):
    """
    Parse SINTAX output and generate summary statistics.
    
    Uses column 3 (confident taxonomy at the specified cutoff) rather than
    column 1 (full taxonomy at any confidence) to avoid inflating counts
    with low-confidence guesses.
    
    SINTAX output format (tab-separated):
    col0: query_id
    col1: full taxonomy with bootstrap values at ALL levels
    col2: strand (+/-)
    col3: taxonomy with only levels >= confidence cutoff
    """
    if not Path(sintax_file).exists():
        return None
    
    stats = {
        'total_otus': 0,
        'assigned': 0,
        'unassigned': 0,
        'assigned_beyond_domain': 0,
        'phylum_counts': {},
        'class_counts': {},
        'order_counts': {},
        'family_counts': {},
        'genus_counts': {},
        'species_counts': {}
    }
    
    with open(sintax_file, 'r') as f:
        for line in f:
            if line.strip() == '':
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            
            stats['total_otus'] += 1
            
            # Use column 3 (confident taxonomy >= cutoff) instead of column 1
            confident_taxonomy = parts[3].strip() if len(parts) > 3 else ""
            
            if confident_taxonomy:
                stats['assigned'] += 1
                
                # Parse confident taxonomy levels (format: d:Eukaryota,p:Phylum,...)
                tax_levels = {}
                for item in confident_taxonomy.split(','):
                    if ':' in item:
                        level, name = item.split(':', 1)
                        tax_levels[level] = name.strip()
                
                # Track OTUs resolved beyond just domain
                if any(k in tax_levels for k in ('p', 'c', 'o', 'f', 'g', 's')):
                    stats['assigned_beyond_domain'] += 1
                
                # Count at each level
                if 'p' in tax_levels:
                    phylum = tax_levels['p']
                    stats['phylum_counts'][phylum] = stats['phylum_counts'].get(phylum, 0) + 1
                if 'c' in tax_levels:
                    cls = tax_levels['c']
                    stats['class_counts'][cls] = stats['class_counts'].get(cls, 0) + 1
                if 'o' in tax_levels:
                    order = tax_levels['o']
                    stats['order_counts'][order] = stats['order_counts'].get(order, 0) + 1
                if 'f' in tax_levels:
                    family = tax_levels['f']
                    stats['family_counts'][family] = stats['family_counts'].get(family, 0) + 1
                if 'g' in tax_levels:
                    genus = tax_levels['g']
                    stats['genus_counts'][genus] = stats['genus_counts'].get(genus, 0) + 1
                if 's' in tax_levels:
                    species = tax_levels['s']
                    stats['species_counts'][species] = stats['species_counts'].get(species, 0) + 1
            else:
                stats['unassigned'] += 1
    
    return stats

def print_taxonomy_summary(stats, marker):
    """Print taxonomy assignment summary statistics."""
    if stats is None:
        print(f"\n[{marker}] No taxonomy results available")
        return
    
    total = stats['total_otus']
    print(f"\n{'='*60}")
    print(f"TAXONOMY SUMMARY: {marker} (confident assignments only)")
    print(f"{'='*60}")
    print(f"Total OTUs: {total}")
    print(f"Assigned (any confident level): {stats['assigned']} ({100*stats['assigned']/total:.1f}%)")
    print(f"Resolved beyond domain: {stats.get('assigned_beyond_domain', 0)} ({100*stats.get('assigned_beyond_domain', 0)/total:.1f}%)")
    print(f"Domain only or unassigned: {total - stats.get('assigned_beyond_domain', 0)} ({100*(total - stats.get('assigned_beyond_domain', 0))/total:.1f}%)")
    
    print(f"\nTop 10 Phyla:")
    for phylum, count in sorted(stats['phylum_counts'].items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"  {phylum}: {count}")
    
    print(f"\nTop 10 Families:")
    for family, count in sorted(stats['family_counts'].items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"  {family}: {count}")
    
    print(f"\nTop 10 Genera:")
    for genus, count in sorted(stats['genus_counts'].items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"  {genus}: {count}")

def main():
    parser = argparse.ArgumentParser(
        description="Assign taxonomy to OTU consensus sequences using VSEARCH SINTAX",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process both 18S and COI
  python scripts/5_assign_taxonomy.py \\
      --input_dir out/Water_eDNA_18S_COI_14_01_26 \\
      --db_18S refs/silva_18s_v123.udb \\
      --db_COI refs/midori2_COI.udb \\
      --threads 12
  
  # Process 18S only with SILVA database
  python scripts/5_assign_taxonomy.py \\
      --input_dir out/Water_eDNA_18S_COI_14_01_26 \\
      --db_18S refs/silva_18s_v123.udb \\
      --threads 12

Database Preparation:
  1. Download SILVA 18S database:
     wget https://www.arb-silva.de/fileadmin/silva_databases/release_123/Exports/SILVA_123_SSURef_tax_silva.fasta.gz
     gunzip SILVA_123_SSURef_tax_silva.fasta.gz
  
  2. Convert to VSEARCH UDB format:
     vsearch --makeudb_usearch SILVA_123_SSURef_tax_silva.fasta --output silva_18s_v123.udb
  
  3. For COI, download MIDORI2 and convert similarly
        """)
    parser.add_argument("--input_dir", required=True, 
                        help="Input directory containing temp_clustering folder with consensus sequences")
    parser.add_argument("--db_18S", 
                        help="Path to 18S rRNA reference database (.udb format required)")
    parser.add_argument("--db_COI", 
                        help="Path to COI reference database (.udb format required)")
    parser.add_argument("--db_JEDI",
                        help="Path to JEDI rRNA reference database (.udb format). "
                             "JEDI targets SSU rRNA V4-V5 (all domains), use SILVA or PR2.")
    parser.add_argument("--threads", type=int, default=4, 
                        help="Number of threads for VSEARCH (default: 4)")
    parser.add_argument("--confidence", type=float, default=0.8, 
                        help="SINTAX confidence threshold, 0-1 (default: 0.8)")
    parser.add_argument("--tag", default=None,
                        help="Subfolder name for taxonomy output (e.g. silva). "
                             "Auto-derived from DB filenames if omitted. Lets multiple DB runs coexist.")
    
    args = parser.parse_args()
    
    # Validate that at least one database is provided
    if not args.db_18S and not args.db_COI and not args.db_JEDI:
        parser.error("At least one database (--db_18S, --db_COI, or --db_JEDI) must be provided")
    
    vsearch_path = find_tool("vsearch")
    
    if not vsearch_path:
        print("ERROR: VSEARCH not found in PATH.", file=sys.stderr)
        print("Please install VSEARCH: conda install -c bioconda vsearch", file=sys.stderr)
        sys.exit(1)
    
    input_dir = Path(args.input_dir)
    temp_dir = input_dir / "temp_clustering"
    
    if not temp_dir.exists():
        print(f"ERROR: temp_clustering directory not found: {temp_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Per-marker DB labels for output paths: taxonomy/{MARKER}/{db_label}/
    # Each marker gets its own subfolder so different DB combos don't duplicate files.
    db_labels = {}
    for marker_name, db_path in [("18S", args.db_18S), ("COI", args.db_COI), ("JEDI", args.db_JEDI)]:
        if db_path:
            db_labels[marker_name] = args.tag or label_from_path(db_path)
    print(f"[tag] Per-marker DB labels: {db_labels}")
    
    markers = {
        "18S": args.db_18S,
        "COI": args.db_COI,
        "JEDI": args.db_JEDI,
    }
    
    print(f"{'='*60}")
    print("TAXONOMY ASSIGNMENT WITH VSEARCH SINTAX")
    print(f"{'='*60}")
    
    for marker, db_path in markers.items():
        if db_path is None:
            print(f"\n[SKIP] No database provided for {marker}")
            continue
        
        print(f"\n{'='*60}")
        print(f"PROCESSING MARKER: {marker}")
        print(f"{'='*60}")
        
        # Input: clean consensus sequences (after chimera removal)
        consensus_clean = temp_dir / f"consensus_{marker}_clean.fasta"
        
        if not consensus_clean.exists():
            print(f"WARNING: Clean consensus file not found: {consensus_clean}")
            print(f"Skipping {marker}. Run clustering first.")
            continue
        
        # Output: taxonomy/{MARKER}/{db_label}/taxonomy_{MARKER}.txt
        db_label = db_labels.get(marker, "default")
        taxonomy_dir = input_dir / "taxonomy" / marker / db_label
        taxonomy_dir.mkdir(parents=True, exist_ok=True)
        taxonomy_file = taxonomy_dir / f"taxonomy_{marker}.txt"
        
        # Run SINTAX
        success = run_sintax(
            consensus_clean,
            taxonomy_file,
            db_path,
            vsearch_path,
            args.threads,
            args.confidence
        )
        
        if success:
            # Parse and summarize results
            stats = parse_sintax_output(taxonomy_file, marker)
            print_taxonomy_summary(stats, marker)
    
    print(f"\n{'='*60}")
    print("TAXONOMY ASSIGNMENT COMPLETE")
    print(f"{'='*60}")
    for m, lab in db_labels.items():
        print(f"  taxonomy/{m}/{lab}/taxonomy_{m}.txt")
    print("\nFiles created:")
    if args.db_18S:
        print("  - taxonomy_18S.txt")
    if args.db_COI:
        print("  - taxonomy_COI.txt")
    print("\nDatabases used:")
    if args.db_18S:
        print(f"  18S: {args.db_18S}")
    if args.db_COI:
        print(f"  COI: {args.db_COI}")
    if args.db_JEDI:
        print(f"  JEDI: {args.db_JEDI}")
    print("See README.md for database setup instructions.")

if __name__ == "__main__":
    main()