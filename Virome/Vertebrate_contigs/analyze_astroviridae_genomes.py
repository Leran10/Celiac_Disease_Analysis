#!/usr/bin/env python3
"""
Comprehensive analysis of Astroviridae contigs:
1. Re-run BLAST to get alignment coordinates
2. Download annotated reference genomes
3. Group contigs by reference genome and check for overlaps
"""

import os
import time
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

# Set your email for NCBI Entrez (required by NCBI)
Entrez.email = "user@institute.edu"

def blast_with_coordinates(seq_record, max_hits=1):
    """
    Run BLASTN and extract alignment coordinates
    """
    print(f"Running BLASTN for {seq_record.id} (length: {len(seq_record.seq)} bp)...")

    try:
        # Submit BLAST query
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            sequence=str(seq_record.seq),
            hitlist_size=max_hits,
            megablast=True
        )

        # Parse BLAST results
        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)

        hits = []
        for alignment in blast_record.alignments[:max_hits]:
            hsp = alignment.hsps[0]  # Get the best HSP

            hit_info = {
                'query_id': seq_record.id,
                'query_length': len(seq_record.seq),
                'hit_accession': alignment.accession,
                'hit_description': alignment.hit_def,
                'hit_length': alignment.length,
                'e_value': hsp.expect,
                'identity': hsp.identities,
                'align_length': hsp.align_length,
                'percent_identity': round((hsp.identities / hsp.align_length) * 100, 2),
                'coverage': round((hsp.align_length / len(seq_record.seq)) * 100, 2),
                'query_start': hsp.query_start,
                'query_end': hsp.query_end,
                'subject_start': hsp.sbjct_start,
                'subject_end': hsp.sbjct_end,
                'strand': 'plus' if hsp.sbjct_start < hsp.sbjct_end else 'minus'
            }
            hits.append(hit_info)

        result_handle.close()
        return hits

    except Exception as e:
        print(f"Error processing {seq_record.id}: {e}")
        return []

def download_genbank_record(accession, output_dir="reference_genomes"):
    """
    Download GenBank file for an accession number
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{accession}.gb")

    # Check if already downloaded
    if os.path.exists(output_file):
        print(f"  {accession} already downloaded, skipping...")
        return output_file

    try:
        print(f"  Downloading {accession}...")
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        genbank_data = handle.read()
        handle.close()

        with open(output_file, 'w') as f:
            f.write(genbank_data)

        print(f"  Saved to {output_file}")
        time.sleep(0.5)  # Be nice to NCBI
        return output_file

    except Exception as e:
        print(f"  Error downloading {accession}: {e}")
        return None

def check_contig_overlap(contig1, contig2):
    """
    Check if two contigs overlap on the reference genome
    """
    # Get subject coordinates
    c1_start = min(contig1['subject_start'], contig1['subject_end'])
    c1_end = max(contig1['subject_start'], contig1['subject_end'])
    c2_start = min(contig2['subject_start'], contig2['subject_end'])
    c2_end = max(contig2['subject_start'], contig2['subject_end'])

    # Check for overlap
    overlap = not (c1_end < c2_start or c2_end < c1_start)

    if overlap:
        overlap_start = max(c1_start, c2_start)
        overlap_end = min(c1_end, c2_end)
        overlap_length = overlap_end - overlap_start
        return True, overlap_length
    else:
        return False, 0

def main():
    input_fasta = "Astroviridae_Celiac_contigs.fas"
    output_coords = "Astroviridae_with_coordinates.tsv"

    print("=" * 70)
    print("STEP 1: Re-running BLAST to get alignment coordinates")
    print("=" * 70)

    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    print(f"Found {len(sequences)} sequences\n")

    all_results = []
    for i, seq_record in enumerate(sequences, 1):
        print(f"[{i}/{len(sequences)}] Processing {seq_record.id}...")
        hits = blast_with_coordinates(seq_record, max_hits=1)
        all_results.extend(hits)

        # Be nice to NCBI servers
        if i < len(sequences):
            print("Waiting 10 seconds before next query...\n")
            time.sleep(10)

    # Save results with coordinates
    print(f"\nSaving results to {output_coords}...")
    df = pd.DataFrame(all_results)
    df.to_csv(output_coords, sep='\t', index=False)
    print(f"Done! Saved {len(all_results)} results\n")

    print("=" * 70)
    print("STEP 2: Downloading annotated reference genomes")
    print("=" * 70)

    # Get unique reference accessions
    unique_refs = df['hit_accession'].unique()
    print(f"Found {len(unique_refs)} unique reference genomes\n")

    for accession in unique_refs:
        download_genbank_record(accession)

    print("\n" + "=" * 70)
    print("STEP 3: Grouping contigs by reference genome")
    print("=" * 70)

    # Group by reference genome
    grouped = df.groupby('hit_accession')

    for ref_acc, group in grouped:
        print(f"\nReference: {ref_acc}")
        print(f"  {group.iloc[0]['hit_description']}")
        print(f"  Number of contigs mapping: {len(group)}")

        if len(group) > 1:
            print(f"  Contigs: {', '.join(group['query_id'].tolist())}")

            # Check for overlaps between contigs
            contigs_list = group.to_dict('records')
            print("\n  Checking for overlaps:")
            for i, c1 in enumerate(contigs_list):
                for j, c2 in enumerate(contigs_list[i+1:], i+1):
                    overlap, overlap_len = check_contig_overlap(c1, c2)
                    if overlap:
                        print(f"    ⚠️  {c1['query_id']} and {c2['query_id']} OVERLAP by {overlap_len} bp")
                    else:
                        gap = min(abs(c1['subject_start'] - c2['subject_end']),
                                abs(c2['subject_start'] - c1['subject_end']))
                        print(f"    ✓  {c1['query_id']} and {c2['query_id']} are separate (gap: {gap} bp)")

            # Show alignment positions
            print("\n  Alignment positions on reference:")
            for _, row in group.iterrows():
                print(f"    {row['query_id']:15} -> ref:{row['subject_start']:6}-{row['subject_end']:6} ({row['strand']})")
        else:
            print(f"  Contig: {group.iloc[0]['query_id']}")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total contigs: {len(df)}")
    print(f"Unique reference genomes: {len(unique_refs)}")
    print(f"Results with coordinates: {output_coords}")
    print(f"Reference genomes downloaded to: reference_genomes/")

if __name__ == "__main__":
    main()
