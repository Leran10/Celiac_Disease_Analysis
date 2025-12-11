#!/usr/bin/env python3
"""
Convert all GenBank files to FASTA format (genome sequences only)
"""

import os
import glob
from Bio import SeqIO

def convert_genbank_to_fasta(genbank_file, output_dir="genome_sequences"):
    """
    Convert a GenBank file to FASTA format
    """
    try:
        # Read GenBank file
        record = SeqIO.read(genbank_file, "genbank")

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Determine output filename
        basename = os.path.basename(genbank_file)
        accession = basename.replace('.gb', '').replace('.txt', '')
        output_file = os.path.join(output_dir, f"{accession}_genome.fasta")

        # Write FASTA file
        SeqIO.write(record, output_file, "fasta")

        print(f"✓ {accession:12} -> {output_file}")
        print(f"  Description: {record.description}")
        print(f"  Length: {len(record.seq):,} bp\n")

        return output_file

    except Exception as e:
        print(f"✗ Error processing {genbank_file}: {e}\n")
        return None

def main():
    # Find all GenBank files
    genbank_dir = "Astroviridae_reference_genomes"
    genbank_files = glob.glob(os.path.join(genbank_dir, "*.gb")) + \
                    glob.glob(os.path.join(genbank_dir, "*.gb.txt"))

    print("=" * 70)
    print(f"Converting {len(genbank_files)} GenBank files to FASTA format")
    print("=" * 70)
    print()

    output_dir = os.path.join(genbank_dir, "genome_sequences")

    converted_files = []
    for gb_file in sorted(genbank_files):
        fasta_file = convert_genbank_to_fasta(gb_file, output_dir)
        if fasta_file:
            converted_files.append(fasta_file)

    print("=" * 70)
    print(f"Summary: {len(converted_files)}/{len(genbank_files)} files converted successfully")
    print(f"Output directory: {output_dir}")
    print("=" * 70)

if __name__ == "__main__":
    main()
