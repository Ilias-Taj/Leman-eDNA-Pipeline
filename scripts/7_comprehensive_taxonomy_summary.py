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

Changelog:
    2025-05-12  Reorder: now runs as step 6 (before BLAST).
                Creates CSV with blank NCBI_TopHit/Identity/Evalue columns.
                BLAST (step 7) fills them in via --update_summary.
                Cosmetic whitespace cleanup only in code.
"""

import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import argparse
import sys
import time
from pathlib import Path
import re

# Ensure scripts/ is on the import path so utils.py can be found from any working directory
sys.path.insert(0, str(Path(__file__).resolve().parent))
from db_tag import DB_RANK_MAPS, label_from_path

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

def parse_taxon_name(raw_taxon, rank_name, db_label):
    """
    Process one taxon name.

    Returns:
        taxon_name
        NCBI TaxID, when provided by MIDORI2
        whether the rank was artificially filled by MIDORI2
    """
    taxon_name = raw_taxon.strip().rstrip("_")
    taxid = ""
    is_imputed = False

    if db_label == "midori2":
        # MIDORI2 names can end with an NCBI TaxID:
        # Arthropoda_6656 -> name=Arthropoda, taxid=6656
        taxid_match = re.match(r"^(.*)_([0-9]+)$", taxon_name)

        if taxid_match:
            taxon_name = taxid_match.group(1).rstrip("_")
            taxid = taxid_match.group(2)

        # MIDORI2 fills missing ranks with names such as:
        # class_Crocodylia
        # order_Vannellidae
        placeholder_prefix = f"{rank_name}_"

        if taxon_name.lower().startswith(placeholder_prefix):
            is_imputed = True
            taxon_name = ""
            taxid = ""

    return taxon_name, taxid, is_imputed

def parse_sintax_taxonomy(taxonomy_file, db_label):
    """
    Parse VSEARCH SINTAX output using the rank structure of the
    selected reference database.
    """
    db_label = db_label.lower()

    if db_label not in DB_RANK_MAPS:
        raise ValueError(
            f"No taxonomy rank mapping defined for database: {db_label}"
        )

    level_map = DB_RANK_MAPS[db_label]
    taxonomy_dict = {}

    with open(taxonomy_file, "r") as handle:
        for line in handle:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 2:
                continue

            # Remove centroid= and abundance annotations from the query ID.
            raw_header = parts[0].split(";")[0]
            raw_header = raw_header.replace("centroid=", "")
            centroid_id = raw_header.split("|")[0].strip()

            full_taxonomy = parts[1].strip()

            tax_dict = {}

            # Initialise only the ranks belonging to this database.
            for rank_name in dict.fromkeys(level_map.values()):
                tax_dict[rank_name] = ""
                tax_dict[f"{rank_name}_conf"] = ""
                tax_dict[f"{rank_name}_taxid"] = ""
                tax_dict[f"{rank_name}_imputed"] = False

            # First split the SINTAX taxonomy into individual ranks.
            for item in full_taxonomy.split(","):
                item = item.strip()

                # Examples:
                # p:Arthropoda(0.95)
                # k:Metazoa_(Animalia)(0.98)
                match = re.match(
                    r"^([dkpcofgs]):(.+)\(([0-9]*\.?[0-9]+)\)$",
                    item,
                )

                if not match:
                    continue

                level_code = match.group(1)
                raw_taxon = match.group(2).strip()
                confidence = float(match.group(3))

                rank_name = level_map.get(level_code)

                if rank_name is None:
                    continue

                taxon_name, taxid, is_imputed = parse_taxon_name(
                    raw_taxon,
                    rank_name,
                    db_label,
                )

                tax_dict[f"{rank_name}_imputed"] = is_imputed

                # Do not report artificial MIDORI2 ranks as assignments.
                if is_imputed:
                    continue

                tax_dict[rank_name] = taxon_name
                tax_dict[f"{rank_name}_conf"] = round(confidence, 2)
                tax_dict[f"{rank_name}_taxid"] = taxid

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

def parse_blast_results_dir(blast_dir, marker):
    """Parse existing BLAST result .txt files from blast_dir for a given marker.

    Reads blast_top*_{abundance,confidence}_{marker}.txt files and returns
    a dict: OTU_ID -> {species, identity, evalue}.
    """
    blast_dir = Path(blast_dir)
    if not blast_dir.is_dir():
        return {}

    results = {}
    # Match files like blast_top10_abundance_COI.txt, blast_top10_confidence_18S.txt
    patterns = [
        f"blast_top*_abundance_{marker}.txt",
        f"blast_top*_confidence_{marker}.txt",
        f"blast_top*_{marker}.txt",  # legacy format without criterion
    ]
    files_found = []
    for pat in patterns:
        files_found.extend(blast_dir.glob(pat))

    for fpath in sorted(set(files_found)):
        with open(fpath, 'r') as f:
            lines = f.readlines()
        reading = False
        for line in lines:
            if line.startswith('---'):
                reading = True
                continue
            if not reading or not line.strip():
                continue
            parts = line.split('|')
            if len(parts) >= 4:
                otu_id = parts[0].strip()
                species = parts[2].strip()
                identity_str = parts[3].strip().replace('%', '')
                evalue_str = parts[4].strip() if len(parts) >= 5 else 'N/A'
                try:
                    identity = float(identity_str) if identity_str and identity_str != '-' else 0
                    results[otu_id] = {
                        'species': species,
                        'identity': round(identity, 1),
                        'evalue': evalue_str if evalue_str else 'N/A',
                    }
                except (ValueError, TypeError):
                    continue
    return results


def detect_db_prefix(db_path, marker):
    """Detect taxonomy column prefix from database filename.

    Inspects the database filename to determine the CSV column prefix:
      silva_18S.udb -> 'SILVA'
      eKOI_COI.udb       -> 'eKOI'
      midori2_COI.udb    -> 'MIDORI2'
    Falls back to 'SILVA' for 18S/JEDI (rRNA), 'eKOI' for COI if no path provided.
    """
    if db_path:
        name = Path(db_path).stem.lower()
        if 'silva' in name:
            return "SILVA"
        elif 'midori' in name:
            return "MIDORI2"
        elif 'ekoi' in name:
            return "eKOI"
        elif 'pr2' in name:
            return "PR2"
        elif 'porter' in name:
            return "Porter"
    # Defaults when no path given
    if marker == "18S":
        return "SILVA"
    return "eKOI"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True, help="Input directory")
    parser.add_argument("--blast_n", type=int, default=0, help="Number of top OTUs to BLAST per marker (0=skip)")
    parser.add_argument("--skip_blast", action='store_true', help="Skip BLAST, use only local taxonomy")
    parser.add_argument("--blast_dir", default=None,
                        help="Path to existing BLAST results directory. When --skip_blast is set, "
                             "existing results from this directory are incorporated into the summary.")
    parser.add_argument("--markers", default=None,
                        help="Comma-separated list of markers (default: auto-detect from merged/ files)")
    # Note: confidence filtering is now done in analysis notebooks, not here
    parser.add_argument("--db_18S", default=None, help="Path to 18S database (for prefix detection)")
    parser.add_argument("--db_COI", default=None, help="Path to COI database (for prefix detection)")
    parser.add_argument("--db_JEDI", default=None, help="Path to JEDI database (for prefix detection)")
    parser.add_argument("--tag", default=None,
                        help="Subfolder name under taxonomy/ and taxonomy_summary/ "
                             "(e.g. 'silva'). Auto-derived from DB filenames "
                             "if omitted. Must match the tag used for script 5.")
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)

    # Per-marker DB labels: taxonomy/{MARKER}/{db_label}/ layout
    db_labels = {}
    for mname, db_path in [("18S", args.db_18S), ("COI", args.db_COI), ("JEDI", args.db_JEDI)]:
        if db_path:
            db_labels[mname] = args.tag or label_from_path(db_path)
    print(f"[tag] Per-marker DB labels: {db_labels}")
    
    # Determine markers to process
    if args.markers:
        markers_to_process = [m.strip().upper() for m in args.markers.split(",")]
    else:
        # Auto-detect from existing abundance files
        merged_dir = input_dir / "merged"
        markers_to_process = []
        for candidate in ["18S", "COI", "JEDI"]:
            if candidate in db_labels and (merged_dir / f"otu_relative_abundance_{candidate}.csv").exists():
                markers_to_process.append(candidate)
        if not markers_to_process:
            print("ERROR: no markers to process (no matching DB + abundance file)", file=sys.stderr)
            sys.exit(1)
    
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
        db_label = db_labels.get(marker)
        if not db_label:
            print(f"  [SKIP] No database provided for {marker}")
            continue
        taxonomy_dir_in = input_dir / "taxonomy" / marker / db_label
        taxonomy_file = taxonomy_dir_in / f"taxonomy_{marker}.txt"
        consensus_file = input_dir / f"temp_clustering/consensus_{marker}_clean.fasta"
        otu_assignment_file = input_dir / f"global_otu_assignment_{marker}.txt"
        
        # Check files exist
        if not abundance_file.exists():
            print(f"  [WARN] Abundance file not found: {abundance_file}")
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
        
        # 3. Load taxonomy (SILVA/PR2 for 18S and JEDI, eKOI for COI)
        # Detect DB prefix from database path
        db_arg = getattr(args, f'db_{marker}', None)
        db_prefix = detect_db_prefix(db_arg, marker)
        print(f"[3/5] Loading local taxonomy assignments ({db_prefix})...")
        taxonomy_assignments = {}

        if taxonomy_file.exists():
            taxonomy_assignments = parse_sintax_taxonomy(
                taxonomy_file,
                db_label,
            )
            print(
                f"  Loaded {len(taxonomy_assignments)} "
                "taxonomy assignments (all confidence levels)"
            )
        else:
            print("  [WARN] Taxonomy file not found")

        
        # 4. Run BLAST if requested, or load existing results from blast_dir
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
        elif args.blast_dir:
            blast_dir = Path(args.blast_dir)
            print(f"[4/5] Loading existing BLAST results from {blast_dir}")
            blast_results = parse_blast_results_dir(blast_dir, marker)
            if blast_results:
                print(f"  Loaded {len(blast_results)} BLAST results for {marker}")
            else:
                print(f"  No existing BLAST results found for {marker} in {blast_dir}")
        else:
            print("[4/5] Skipping BLAST (no --blast_dir provided)")
        
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
            
            # Add taxonomy (SILVA/PR2 for 18S and JEDI, MIDORI/eKOI for COI)
            # All levels stored with confidence — filtering done in notebooks
            # With isONclust3, OTU ID is the SINTAX query ID directly.
            # Fall back to centroid UUID lookup for old VSEARCH-based runs.
            rank_map = DB_RANK_MAPS[db_label]
            tax_levels = [
                rank_name.title()
                for rank_name in dict.fromkeys(
                    DB_RANK_MAPS[db_label].values()
                )
            ]
            
            centroid_id = otu_to_centroid.get(otu_id, '')
            tax_key = otu_id if otu_id in taxonomy_assignments else centroid_id
            if tax_key in taxonomy_assignments:
                tax = taxonomy_assignments[tax_key]
                for level in tax_levels:
                    level_key = level.lower()

                    row[f"{db_prefix}_{level}"] = tax.get(
                        level_key,
                        "",
                    )
                    row[f"{db_prefix}_{level}_Conf"] = tax.get(
                        f"{level_key}_conf",
                        "",
                    )

                    if db_label == "midori2":
                        row[f"{db_prefix}_{level}_TaxID"] = tax.get(
                            f"{level_key}_taxid",
                            "",
                        )
                        row[f"{db_prefix}_{level}_Imputed"] = tax.get(
                            f"{level_key}_imputed",
                            False,
                        )

            else:
                for level in tax_levels:
                    row[f"{db_prefix}_{level}"] = ""
                    row[f"{db_prefix}_{level}_Conf"] = ""

                    if db_label == "midori2":
                        row[f"{db_prefix}_{level}_TaxID"] = ""
                        row[f"{db_prefix}_{level}_Imputed"] = False

            

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
        output_dir = input_dir / "taxonomy_summary" / marker / db_label
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"comprehensive_taxonomy_{marker}.csv"
        summary_df.to_csv(output_file, index=False)
        print(f"\n[OK] Saved: {output_file}")
        print(f"  {len(summary_df)} OTUs with taxonomy and abundance data")
        
        # Print summary statistics
        print(f"\nSummary Statistics for {marker}:")
        print(f"  Total OTUs: {len(summary_df)}")
        
        broad_rank = {
            "pr2": "Division",
            "silva": "Phylum",
            "midori2": "Phylum",
            "porter": "Phylum",
        }[db_label]

        assigned_column = f"{db_prefix}_{broad_rank}"

        assigned = summary_df[
            summary_df[assigned_column].fillna("") != ""
        ].shape[0]
        
        print(f"  With {db_prefix} taxonomy: {assigned} ({100*assigned/len(summary_df):.1f}%)")
        if blast_results:
            blasted = summary_df[summary_df['NCBI_TopHit'] != ''].shape[0]
            print(f"  With NCBI BLAST: {blasted} ({100*blasted/len(summary_df):.1f}%)")
    
    print("\n" + "=" * 80)
    print("COMPREHENSIVE TAXONOMY SUMMARY COMPLETE")
    for m, lab in db_labels.items():
        print(f"  taxonomy_summary/{m}/{lab}/comprehensive_taxonomy_{m}.csv")
    print("=" * 80)

if __name__ == "__main__":
    main()
