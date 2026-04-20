#!/usr/bin/env python3
"""
Convert PR2 SSU DADA2 fasta to VSEARCH SINTAX format for 18S classification.

PR2 DADA2 format (header):
  >Accession Kingdom;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species

VSEARCH SINTAX format:
  >Accession;tax=d:Kingdom,k:Supergroup,p:Division,c:Subdivision,o:Class,f:Family,g:Genus,s:Species;

Mapping (PR2 9-level → SINTAX 8-level):
  d: = Kingdom (Eukaryota)
  k: = Supergroup (Obazoa, TSAR, Archaeplastida...)
  p: = Division (Opisthokonta, Alveolata, Stramenopiles...)
  c: = Subdivision (Metazoa, Ciliophora, Dinoflagellata...)
  o: = Class (used as Order in SINTAX for broader classification)
  f: = Family (PR2 Order+Family collapsed when _X placeholders)
  g: = Genus
  s: = Species

Note: PR2 has 9 taxonomic levels but SINTAX supports 8. We collapse
Order into o: when it's a placeholder, otherwise keep Class as o: and
Order as f: prefix.

Usage:
  python refs/convert_pr2_to_sintax.py \
      --input refs/pr2_version_5.1.1_SSU_dada2.fasta \
      --output refs/pr2_18S_SINTAX.fasta

Then build the .udb:
  ./env/bin/vsearch --makeudb_usearch refs/pr2_18S_SINTAX.fasta \
      --output refs/pr2_18S_v511.udb
"""

import argparse
import sys
import re


def clean_field(field):
    """Remove trailing _X, _XX, _XXX placeholders and clean field."""
    if field is None:
        return ""
    field = field.strip()
    # Remove placeholder suffixes like _X, _XX, _XXX
    field = re.sub(r'_X+$', '', field)
    if not field:
        return ""
    # Replace spaces with underscores for SINTAX compatibility
    field = field.replace(' ', '_')
    return field


def convert_header(header_line):
    """Convert a PR2 DADA2 header to SINTAX format.

    Expected format: >Accession Kingdom;Supergroup;Division;...;Species
    """
    line = header_line.lstrip('>').strip()

    # Split on first space: accession vs taxonomy
    parts = line.split(' ', 1)
    if len(parts) < 2:
        return None

    accession = parts[0].strip()
    taxonomy = parts[1].strip()

    # Split taxonomy by semicolons
    fields = [f.strip() for f in taxonomy.split(';')]

    # PR2 has 9 levels: Kingdom;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species
    if len(fields) < 7:
        return None

    # Pad to 9 fields if shorter
    while len(fields) < 9:
        fields.append('')

    kingdom     = clean_field(fields[0])
    supergroup  = clean_field(fields[1])
    division    = clean_field(fields[2])
    subdivision = clean_field(fields[3])
    pr2_class   = clean_field(fields[4])
    pr2_order   = clean_field(fields[5])
    pr2_family  = clean_field(fields[6])
    genus       = clean_field(fields[7]) if len(fields) > 7 else ''
    species     = clean_field(fields[8]) if len(fields) > 8 else ''

    # Build SINTAX taxonomy string (8 levels max)
    # Use class as order-level (o:) and family as family-level (f:)
    # If order is empty/placeholder, use class as order
    order_name = pr2_order if pr2_order else pr2_class
    family_name = pr2_family if pr2_family else pr2_order

    tax_parts = []
    if kingdom:
        tax_parts.append(f"d:{kingdom}")
    if supergroup:
        tax_parts.append(f"k:{supergroup}")
    if division:
        tax_parts.append(f"p:{division}")
    if subdivision:
        tax_parts.append(f"c:{subdivision}")
    if order_name:
        tax_parts.append(f"o:{order_name}")
    if family_name:
        tax_parts.append(f"f:{family_name}")
    if genus:
        tax_parts.append(f"g:{genus}")
    if species:
        tax_parts.append(f"s:{species}")

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
        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                new_header = convert_header(line)
                if new_header:
                    fout.write(new_header + '\n')
                    converted += 1
                else:
                    skipped += 1
            else:
                # Only write sequence if the previous header was valid
                if converted > 0 or skipped == 0:
                    fout.write(line + '\n')

    print(f"Converted {converted:,} sequences ({skipped:,} skipped)")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
