#!/usr/bin/env python3
"""
Extract genome sequence from GenBank file
"""

from Bio import SeqIO

def extract_sequence(genbank_file):
    """
    Extract the genome sequence from a GenBank file
    """
    record = SeqIO.read(genbank_file, "genbank")

    print(f"Accession: {record.id}")
    print(f"Description: {record.description}")
    print(f"Genome length: {len(record.seq)} bp")
    print(f"\nFirst 200 bp of genome:")
    print(record.seq[:200])
    print(f"\nLast 200 bp of genome:")
    print(record.seq[-200:])
    print(f"\nFull sequence (first 500 characters):")
    print(str(record.seq)[:500])

    # Save full sequence to FASTA
    output_file = genbank_file.replace('.gb', '_genome.fasta').replace('.txt', '')
    SeqIO.write(record, output_file, "fasta")
    print(f"\nFull genome saved to: {output_file}")

    return record

if __name__ == "__main__":
    genbank_file = "reference_genomes/KF039910.gb.txt"
    record = extract_sequence(genbank_file)
