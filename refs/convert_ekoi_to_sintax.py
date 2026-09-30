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
import re


def clean_taxon(value, rank):
    """Clean a taxon name and omit unresolved placeholder ranks."""
    value = value.strip()

    if not value:
        return ""

    # eKOI/PR2 placeholders such as Tardigrada_X or Tardigrada_XX
    # indicate an unresolved rank, not a real order or family.
    if re.search(r"_X+$", value):
        return ""

    # Names such as Milnesium_sp. are not species-level assignments.
    if rank == "s" and re.search(r"(?:^|[_ ])sp\.?$", value, re.IGNORECASE):
        return ""

    # VSEARCH does not permit commas or semicolons in taxon names.
    value = value.replace(" ", "_")
    value = value.replace(",", "_").replace(";", "_")

    return value


def convert_header(header_line):
    """Convert an eKOI PR2-format header to VSEARCH SINTAX format."""
    line = header_line.lstrip(">").strip().rstrip(";")
    fields = [field.strip() for field in line.split(";")]

    if len(fields) != 10:
        print(
            f"WARNING: Expected 10 fields, found {len(fields)}: "
            f"{header_line.rstrip()}",
            file=sys.stderr,
        )
        return None

    accession = fields[9]
    if not accession:
        print(
            f"WARNING: Header has no accession: {header_line.rstrip()}",
            file=sys.stderr,
        )
        return None

    # eKOI PR2 hierarchy:
    # Domain;Supergroup;Division;Subdivision;Class;Order;
    # Family;Genus;Species;Accession
    #
    # Subdivision (fields[3]) is deliberately omitted because SINTAX
    # provides only eight rank codes.
    rank_values = (
        ("d", clean_taxon(fields[0], "d")),  # Domain
        ("k", clean_taxon(fields[1], "k")),  # Supergroup
        ("p", clean_taxon(fields[2], "p")),  # Division
        ("c", clean_taxon(fields[4], "c")),  # Class
        ("o", clean_taxon(fields[5], "o")),  # Order
        ("f", clean_taxon(fields[6], "f")),  # Family
        ("g", clean_taxon(fields[7], "g")),  # Genus
        ("s", clean_taxon(fields[8], "s")),  # Species
    )

    tax_parts = [
        f"{rank}:{value}"
        for rank, value in rank_values
        if value
    ]

    if not tax_parts:
        print(
            f"WARNING: Header has no usable taxonomy: {header_line.rstrip()}",
            file=sys.stderr,
        )
        return None

    accession = accession.replace(" ", "_").replace(";", "_")
    return f">{accession};tax={','.join(tax_parts)};"


def main():
    parser = argparse.ArgumentParser(description="Convert eKOI PR2 fasta to VSEARCH SINTAX format")
    parser.add_argument("--input", required=True, help="Input eKOI PR2 fasta file")
    parser.add_argument("--output", required=True, help="Output SINTAX-formatted fasta file")
    args = parser.parse_args()
    
    converted = 0
    skipped = 0
    write_sequence = False

    with (
        open(args.input, "r", encoding="utf-8") as fin,
        open(args.output, "w", encoding="utf-8") as fout,
    ):
        for line in fin:
            if line.startswith(">"):
                new_header = convert_header(line)

                if new_header is None:
                    skipped += 1
                    write_sequence = False
                    continue

                fout.write(new_header + "\n")
                converted += 1
                write_sequence = True

            elif write_sequence:
                fout.write(line)
    
    print(f"Converted {converted} sequences ({skipped} skipped)")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
