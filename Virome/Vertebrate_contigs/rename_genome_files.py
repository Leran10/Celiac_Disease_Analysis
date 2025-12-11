#!/usr/bin/env python3
"""
Script to rename reference genome FASTA files with their corresponding contig IDs
"""

import os
import pandas as pd
from collections import defaultdict

def process_family(family_name, base_dir):
    """Process a viral family directory and rename genome files"""

    family_dir = os.path.join(base_dir, family_name)
    coords_file = os.path.join(family_dir, f"{family_name}_with_coordinates.tsv")
    genome_dir = os.path.join(family_dir, f"{family_name}_reference_genomes", "genome_sequences")

    if not os.path.exists(coords_file):
        print(f"Warning: {coords_file} not found, skipping {family_name}")
        return

    if not os.path.exists(genome_dir):
        print(f"Warning: {genome_dir} not found, skipping {family_name}")
        return

    # Read the coordinates file
    df = pd.read_csv(coords_file, sep='\t')

    # Create mapping: accession -> list of contig IDs
    accession_to_contigs = defaultdict(list)
    for _, row in df.iterrows():
        contig_id = row['query_id']
        accession = row['hit_accession']
        accession_to_contigs[accession].append(contig_id)

    # Get list of genome files
    genome_files = [f for f in os.listdir(genome_dir) if f.endswith(('.fasta', '.fasta.txt'))]

    print(f"\nProcessing {family_name}:")
    print(f"  Found {len(genome_files)} genome files")
    print(f"  Found {len(accession_to_contigs)} unique accessions in coordinates file")

    renamed_count = 0
    for genome_file in genome_files:
        # Extract accession from filename (e.g., "PP901925_genome.fasta" -> "PP901925")
        accession = genome_file.split('_genome')[0]

        if accession in accession_to_contigs:
            contigs = accession_to_contigs[accession]
            # Sort contigs to ensure consistent naming
            contigs_sorted = sorted(contigs, key=lambda x: int(x.split('_')[1]))

            # Create new filename with contig IDs
            contig_prefix = "_".join(contigs_sorted)

            # Preserve original file extension
            if genome_file.endswith('.fasta.txt'):
                new_filename = f"{contig_prefix}_{accession}_genome.fasta.txt"
            else:
                new_filename = f"{contig_prefix}_{accession}_genome.fasta"

            old_path = os.path.join(genome_dir, genome_file)
            new_path = os.path.join(genome_dir, new_filename)

            # Rename the file
            os.rename(old_path, new_path)
            print(f"  ✓ {genome_file} -> {new_filename}")
            renamed_count += 1
        else:
            print(f"  ⚠ No contig mapping found for {accession} (file: {genome_file})")

    print(f"  Renamed {renamed_count}/{len(genome_files)} files")

    return renamed_count

def main():
    base_dir = "/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs"

    families = [
        "Adenoviridae",
        "Anelloviridae",
        "Caliciviridae",
        "Picornaviridae"
    ]

    total_renamed = 0
    for family in families:
        renamed = process_family(family, base_dir)
        if renamed:
            total_renamed += renamed

    print(f"\n{'='*60}")
    print(f"Total files renamed across all families: {total_renamed}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
