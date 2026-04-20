#!/usr/bin/env python3
"""
Convert eKOI PR2 fasta to VSEARCH SINTAX format.

eKOI PR2 format (10 semicolon-separated fields):
  >Kingdom;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species;Accession

VSEARCH SINTAX format:
  >Accession;tax=d:Kingdom,k:Supergroup,p:Division,c:Subdivision,o:Class_Order,f:Family,g:Genus,s:Species;

Mapping (PR2 9-level → SINTAX 8-level):
  d: = Kingdom (Eukaryota)
  k: = Supergroup (Obazoa, TSAR, Archaeplastida...)
  p: = Division (Opisthokonta, Alveolata...)  
  c: = Subdivision (Metazoa, Fungi, Apicomplexa...)
  o: = Class (Tardigrada, Nemertea...)
  f: = Family (Order+Family collapsed when _X placeholders)
  g: = Genus
  s: = Species

Usage:
  python refs/convert_ekoi_to_sintax.py \
      --input refs/eKOI_taxonomy_PR2_ver1.fasta \
      --output refs/eKOI_COI_SINTAX.fasta
"""

import argparse
import sys


def convert_header(header_line):
    """Convert a single eKOI PR2 header to SINTAX format."""
    # Remove leading '>'
    fields = header_line.lstrip('>').strip().split(';')
    
    if len(fields) != 10:
        print(f"WARNING: Unexpected field count ({len(fields)}): {header_line}", file=sys.stderr)
        return None
    
    kingdom    = fields[0]  # Eukaryota
    supergroup = fields[1]  # Obazoa, TSAR, Archaeplastida...
    division   = fields[2]  # Opisthokonta, Alveolata, Rhodophyta...
    subdivision= fields[3]  # Metazoa, Fungi, Apicomplexa...
    pr2_class  = fields[4]  # e.g. Tardigrada, Nemertea...
    pr2_order  = fields[5]  # e.g. Tardigrada_X or real order
    pr2_family = fields[6]  # e.g. Tardigrada_XX or real family
    genus      = fields[7]  # e.g. Milnesium
    species    = fields[8]  # e.g. Milnesium_sp.
    accession  = fields[9]  # GenBank accession
    
    # For SINTAX, resolve placeholder _X, _XX suffixes:
    # If order is just "<Class>_X", use the class name as order
    order_name = pr2_order if not pr2_order.endswith('_X') else pr2_class
    # If family is just "<something>_XX", use a combined class name
    family_name = pr2_family if not pr2_family.endswith('_XX') else pr2_class
    
    # Build SINTAX taxonomy string
    tax = (
        f"d:{kingdom},"
        f"k:{supergroup},"
        f"p:{division},"
        f"c:{subdivision},"
        f"o:{order_name},"
        f"f:{family_name},"
        f"g:{genus},"
        f"s:{species}"
    )
    
    return f">{accession};tax={tax};"


def main():
    parser = argparse.ArgumentParser(description="Convert eKOI PR2 fasta to VSEARCH SINTAX format")
    parser.add_argument("--input", required=True, help="Input eKOI PR2 fasta file")
    parser.add_argument("--output", required=True, help="Output SINTAX-formatted fasta file")
    args = parser.parse_args()
    
    converted = 0
    skipped = 0
    
    with open(args.input, 'r') as fin, open(args.output, 'w') as fout:
        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                new_header = convert_header(line)
                if new_header:
                    fout.write(new_header + '\n')
                    converted += 1
                else:
                    # Write original as fallback
                    fout.write(line + '\n')
                    skipped += 1
            else:
                fout.write(line + '\n')
    
    print(f"Converted {converted} sequences ({skipped} skipped)")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
