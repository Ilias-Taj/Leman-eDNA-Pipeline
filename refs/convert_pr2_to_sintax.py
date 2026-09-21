#!/usr/bin/env python3
"""
Convert PR2 SSU DADA2 fasta to VSEARCH SINTAX format for 18S classification.

PR2 DADA2 format (header) - NO accession, just taxonomy:
  >Kingdom;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species;

VSEARCH SINTAX format:
  >Accession;tax=d:Domain,k:Supergroup,p:Division,c:Class,o:Order,f:Family,g:Genus,s:Species;
  
Mapping (PR2 9-level -> SINTAX 8-level):
  d: = Kingdom (Eukaryota)
  k: = Supergroup (Obazoa, TSAR, Archaeplastida...)
  p: = Division (Opisthokonta, Alveolata, Stramenopiles...)
  c: = Class
  o: = Order
  f: = Family (or Order if Family is placeholder)
  g: = Genus
  s: = Species
  
PR2's Subdivision rank is intentionally omitted because SINTAX provides only eight rank codes.

Usage:
  python3 refs/convert_pr2_to_sintax.py \
      --input refs/pr2_version_5.1.1_SSU_dada2.fasta.gz \
      --output refs/pr2_18S_SINTAX.fasta

Then build the .udb:
  ./env/bin/vsearch --makeudb_usearch refs/pr2_18S_SINTAX.fasta \
      --output refs/pr2_18S_v511.udb
"""

import argparse
import re


def clean_field(field):
    """Remove trailing _X, _XX, _XXX placeholders and clean field."""
    if field is None:
        return ""
    field = field.strip()
    # Remove placeholder suffixes like _X, _XX, _XXX
    field = re.sub(r'_X+$', '', field)
    # Remove _sp. suffix
    field = re.sub(r'_sp\.$', '', field)
    if not field:
        return ""
    # Replace spaces with underscores for SINTAX compatibility
    field = field.replace(' ', '_')
    return field


def convert_header(header_line, seq_idx):
    """Convert a PR2 DADA2 header to SINTAX format.

    PR2 DADA2 format has NO accession - header is just semicolon-separated taxonomy:
      >Kingdom;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species;
    We generate a synthetic accession from the sequence index.
    """
    line = header_line.lstrip('>').strip().rstrip(';')

    # Check if there's an accession (space-separated) or pure taxonomy
    if ' ' in line:
        parts = line.split(' ', 1)
        accession = parts[0].strip()
        taxonomy = parts[1].strip().rstrip(';')
    else:
        # Pure taxonomy (PR2 DADA2 format) - generate accession
        accession = f'PR2_{seq_idx:06d}'
        taxonomy = line

    # Split taxonomy by semicolons
    fields = [f.strip() for f in taxonomy.split(';')]

    # PR2 has 9 levels: Kingdom;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species
    if len(fields) < 7:
        return None

    # Pad to 9 fields if shorter
    while len(fields) < 9:
        fields.append('')

    kingdom      = clean_field(fields[0])
    supergroup  = clean_field(fields[1])
    division    = clean_field(fields[2])
    
    # fields[3] is PR2 Subdivision. SINTAX has no corresponding
    # rank code, so it is deliberately omitted.
    class_name  = clean_field(fields[4])
    order_name  = clean_field(fields[5])
    family_name = clean_field(fields[6])
    genus       = clean_field(fields[7])
    species     = clean_field(fields[8])
    
    rank_values = (
        ("d", kingdom),
        ("k", supergroup),
        ("p", division),
        ("c", class_name),
        ("o", order_name),
        ("f", family_name),
        ("g", genus),
        ("s", species),
    )
    
    tax_parts = [
        f"{rank}:{value}"
        for rank, value in rank_values
        if value
    ]

    if not tax_parts:
        return None

    tax_str = ','.join(tax_parts)
    return f">{accession};tax={tax_str};"


def main():
    parser = argparse.ArgumentParser(
        description="Convert PR2 SSU DADA2 fasta to VSEARCH SINTAX format"
    )
    parser.add_argument("--input", required=True,
                        help="Input PR2 DADA2 fasta (plain or .gz)")
    parser.add_argument("--output", required=True,
                        help="Output SINTAX-formatted fasta file")
    args = parser.parse_args()

    import gzip

    # Support gzipped input
    if args.input.endswith('.gz'):
        fin = gzip.open(args.input, 'rt', encoding='utf-8')
    else:
        fin = open(args.input, 'r', encoding='utf-8')

    converted = 0
    skipped = 0

    with fin, open(args.output, 'w', encoding='utf-8') as fout:
        last_valid = False
        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                new_header = convert_header(line, converted + skipped + 1)
                if new_header:
                    fout.write(new_header + '\n')
                    converted += 1
                    last_valid = True
                else:
                    skipped += 1
                    last_valid = False
            else:
                # Only write sequence if the previous header was valid
                if last_valid:
                    fout.write(line + '\n')

    print(f"Converted {converted:,} sequences ({skipped:,} skipped)")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
