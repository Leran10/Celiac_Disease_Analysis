#!/usr/bin/env python3
"""
Script to add Species and mapping_completeness columns to TSV files.
Extracts species information and genome completeness from hit_description column.
"""

import pandas as pd
import re
from pathlib import Path

def extract_species(hit_description):
    """
    Extract species-level annotation from hit_description.
    Examples:
    - "Human astrovirus 1 isolate..." -> "Human astrovirus 1"
    - "Coxsackievirus A4 strain..." -> "Coxsackievirus A4"
    - "Betatorquevirus sp., isolate..." -> "Betatorquevirus sp."
    """
    if pd.isna(hit_description):
        return None

    # Remove "MAG: " or "MAG TPA_asm: " prefix if present
    desc = re.sub(r'^MAG(\s+TPA_asm)?:\s*', '', hit_description)

    # Common patterns to extract species
    # Pattern 1: Genus species number (e.g., "Human astrovirus 1")
    # Pattern 2: Genus species letter+number (e.g., "Coxsackievirus A4", "Rhinovirus A12")
    # Pattern 3: Genus species word (e.g., "Torque teno virus", "Chicken anemia virus")
    # Pattern 4: Genus sp. (e.g., "Betatorquevirus sp.")
    # Pattern 5: Multiple word genus + species (e.g., "Torque teno mini virus")

    # Try to extract up to the first occurrence of isolate/strain/gene/genomic/polyprotein/putative
    match = re.match(r'^([^,]+?)\s+(?:isolate|strain|gene|genomic|polyprotein|putative|ORF)', desc)
    if match:
        species = match.group(1).strip()
        return species

    # If no match with those keywords, try to extract the first part before comma
    match = re.match(r'^([^,]+)', desc)
    if match:
        species = match.group(1).strip()
        # Remove trailing words like "genome", "sequence", etc.
        species = re.sub(r'\s+(genome|sequence)$', '', species)
        return species

    return hit_description

def extract_completeness(hit_description):
    """
    Extract genome/gene completeness information from hit_description.
    Returns one of: "complete genome", "partial genome", "complete cds",
    "partial cds", "complete ORF", "genomic sequence", or None
    """
    if pd.isna(hit_description):
        return None

    # Convert to lowercase for case-insensitive matching
    desc_lower = hit_description.lower()

    # Check for completeness indicators in order of specificity
    if 'complete genome' in desc_lower:
        return 'complete genome'
    elif 'partial genome' in desc_lower:
        return 'partial genome'
    elif 'complete cds' in desc_lower:
        return 'complete cds'
    elif 'partial cds' in desc_lower:
        return 'partial cds'
    elif 'complete orf' in desc_lower:
        return 'complete ORF'
    elif 'genomic sequence' in desc_lower:
        return 'genomic sequence'
    else:
        return None

def process_tsv_file(file_path):
    """
    Process a single TSV file and add Species and mapping_completeness columns.
    """
    print(f"Processing {file_path}...")

    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')

    # Check if hit_description column exists
    if 'hit_description' not in df.columns:
        print(f"  Warning: 'hit_description' column not found in {file_path}")
        return

    # Extract Species and mapping_completeness
    df['Species'] = df['hit_description'].apply(extract_species)
    df['mapping_completeness'] = df['hit_description'].apply(extract_completeness)

    # Save back to the same file
    df.to_csv(file_path, sep='\t', index=False)

    print(f"  Added Species and mapping_completeness columns")
    print(f"  Total rows: {len(df)}")
    print(f"  Species extracted: {df['Species'].notna().sum()}")
    print(f"  Completeness extracted: {df['mapping_completeness'].notna().sum()}")
    print()

def main():
    """
    Main function to process all TSV files with coordinates.
    """
    # Define the base directory
    base_dir = Path("/Users/leranwang/Handley Lab Dropbox/16S/Celiac/Vertebrate_contigs")

    # Define the virus families
    families = ['Astroviridae', 'Anelloviridae', 'Adenoviridae', 'Caliciviridae', 'Picornaviridae']

    # Process each family
    for family in families:
        tsv_file = base_dir / family / f"{family}_with_coordinates.tsv"

        if tsv_file.exists():
            process_tsv_file(tsv_file)
        else:
            print(f"Warning: {tsv_file} not found, skipping...")

    print("All files processed successfully!")

if __name__ == "__main__":
    main()
