#!/usr/bin/env python3
"""
Extract top 1 hit per contig from BLAST results TSV
"""

import pandas as pd

def extract_top_hits(input_file, output_file, top_n=1):
    """
    Extract top N hits per contig from BLAST results
    """
    # Read the TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Group by query_id and take top N hits
    top_hits = df.groupby('query_id').head(top_n).reset_index(drop=True)

    # Save to new file
    top_hits.to_csv(output_file, sep='\t', index=False)

    print(f"Input file: {input_file}")
    print(f"Total hits: {len(df)}")
    print(f"Unique contigs: {df['query_id'].nunique()}")
    print(f"Top {top_n} hits per contig: {len(top_hits)}")
    print(f"Output saved to: {output_file}")

    return top_hits

if __name__ == "__main__":
    # Extract top 1 hit
    input_file = "Astroviridae_Celiac_contigs.tsv"

    print("=" * 60)
    print("Extracting top 1 hit per contig...")
    print("=" * 60)
    top1_df = extract_top_hits(input_file, "Astroviridae_top1_hits.tsv", top_n=1)

    print("\n" + "=" * 60)
    print("Top 1 hits summary:")
    print("=" * 60)
    print(top1_df[['query_id', 'hit_accession', 'hit_description', 'percent_identity', 'coverage']].to_string(index=False))
