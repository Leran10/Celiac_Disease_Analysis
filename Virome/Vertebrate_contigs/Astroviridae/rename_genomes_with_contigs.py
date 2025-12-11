#!/usr/bin/env python3
"""
Rename genome sequence files to include contig IDs based on BLAST hits
"""

import os
import sys
import pandas as pd
from pathlib import Path
from collections import defaultdict

def rename_genomes_for_family(family_name, base_dir="."):
    """
    Rename genome sequences for a viral family to include contig IDs

    Args:
        family_name: Name of viral family (e.g., 'Adenoviridae')
        base_dir: Base directory containing all family folders
    """
    print(f"\n{'='*70}")
    print(f"Processing {family_name}")
    print(f"{'='*70}")

    # Set paths
    family_dir = Path(base_dir) / family_name
    tsv_file = family_dir / f"{family_name}_with_coordinates.tsv"
    genome_dir = family_dir / f"{family_name}_reference_genomes" / "genome_sequences"

    # Check if files exist
    if not tsv_file.exists():
        print(f"✗ TSV file not found: {tsv_file}")
        return False

    if not genome_dir.exists():
        print(f"✗ Genome directory not found: {genome_dir}")
        return False

    # Read TSV file
    df = pd.read_csv(tsv_file, sep='\t')
    print(f"✓ Loaded {len(df)} contig-to-reference mappings")

    # Group contigs by reference accession
    contig_map = defaultdict(list)
    for _, row in df.iterrows():
        accession = row['hit_accession']
        contig_id = row['query_id']
        contig_map[accession].append(contig_id)

    print(f"✓ Found {len(contig_map)} unique reference genomes")

    # Rename genome files
    renamed_count = 0
    skipped_count = 0

    for accession, contig_ids in sorted(contig_map.items()):
        # Find the genome file (could be .fasta or .fa)
        old_file = genome_dir / f"{accession}_genome.fasta"
        if not old_file.exists():
            old_file = genome_dir / f"{accession}_genome.fa"

        if not old_file.exists():
            print(f"  ⚠ Genome file not found for {accession}")
            skipped_count += 1
            continue

        # Create new filename with contig IDs
        contig_prefix = "_".join(sorted(contig_ids))
        new_file = genome_dir / f"{contig_prefix}_{accession}_genome.fasta"

        # Check if already renamed
        if old_file.name.startswith("contig_"):
            print(f"  ⊙ Already renamed: {old_file.name}")
            continue

        # Rename the file
        try:
            old_file.rename(new_file)
            contig_str = ", ".join(sorted(contig_ids))
            print(f"  ✓ {accession:15} -> {contig_prefix}_{accession}_genome.fasta")
            print(f"    Contigs: {contig_str}")
            renamed_count += 1
        except Exception as e:
            print(f"  ✗ Error renaming {accession}: {e}")
            skipped_count += 1

    print(f"\n{'='*70}")
    print(f"Summary for {family_name}:")
    print(f"  Renamed: {renamed_count} files")
    if skipped_count > 0:
        print(f"  Skipped: {skipped_count} files")
    print(f"{'='*70}")

    return True

def main():
    # Get the base directory (parent of Astroviridae)
    script_dir = Path(__file__).parent
    base_dir = script_dir.parent

    # Viral families to process
    families = ['Adenoviridae', 'Anelloviridae', 'Caliciviridae', 'Picornaviridae']

    if len(sys.argv) > 1:
        # Process specific family if provided as argument
        families = [sys.argv[1]]

    print("\n" + "="*70)
    print("Renaming genome sequences with contig IDs")
    print("="*70)
    print(f"Base directory: {base_dir}")
    print(f"Families to process: {', '.join(families)}")

    success_count = 0
    for family in families:
        if rename_genomes_for_family(family, base_dir):
            success_count += 1

    print("\n" + "="*70)
    print(f"FINAL SUMMARY: {success_count}/{len(families)} families processed successfully")
    print("="*70)

if __name__ == "__main__":
    main()
