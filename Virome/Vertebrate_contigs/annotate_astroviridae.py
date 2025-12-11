#!/usr/bin/env python3
"""
Annotate Astroviridae contigs using NCBI BLASTN
"""

import os
import time
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

# Set your email for NCBI Entrez (required by NCBI)
Entrez.email = "user@institute.edu"

def blast_sequence(seq_record, max_hits=5):
    """
    Run BLASTN for a sequence and return top hits
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
            hsp = alignment.hsps[0]  # Get the best HSP (High-scoring Segment Pair)

            hit_info = {
                'query_id': seq_record.id,
                'query_length': len(seq_record.seq),
                'hit_id': alignment.hit_id,
                'hit_def': alignment.hit_def,
                'hit_accession': alignment.accession,
                'length': alignment.length,
                'e_value': hsp.expect,
                'identity': hsp.identities,
                'align_length': hsp.align_length,
                'percent_identity': round((hsp.identities / hsp.align_length) * 100, 2),
                'coverage': round((hsp.align_length / len(seq_record.seq)) * 100, 2)
            }
            hits.append(hit_info)

        result_handle.close()
        return hits

    except Exception as e:
        print(f"Error processing {seq_record.id}: {e}")
        return []

def main():
    input_file = "Astroviridae_Celiac_contigs.fas"
    output_file = "Astroviridae_Celiac_contigs.tsv"

    print(f"Reading sequences from {input_file}...")
    sequences = list(SeqIO.parse(input_file, "fasta"))
    print(f"Found {len(sequences)} sequences")

    all_results = []

    for i, seq_record in enumerate(sequences, 1):
        print(f"\n[{i}/{len(sequences)}] Processing {seq_record.id}...")
        hits = blast_sequence(seq_record, max_hits=5)
        all_results.extend(hits)

        # Be nice to NCBI servers - wait between queries
        if i < len(sequences):
            print("Waiting 10 seconds before next query...")
            time.sleep(10)

    # Write results to TSV
    print(f"\nWriting results to {output_file}...")
    with open(output_file, 'w') as f:
        # Header
        f.write("query_id\tquery_length\thit_accession\thit_description\thit_length\t"
                "e_value\tidentity\talign_length\tpercent_identity\tcoverage\n")

        # Data
        for hit in all_results:
            f.write(f"{hit['query_id']}\t{hit['query_length']}\t{hit['hit_accession']}\t"
                   f"{hit['hit_def']}\t{hit['length']}\t{hit['e_value']}\t"
                   f"{hit['identity']}\t{hit['align_length']}\t{hit['percent_identity']}\t"
                   f"{hit['coverage']}\n")

    print(f"\nDone! Results written to {output_file}")
    print(f"Total hits: {len(all_results)}")

if __name__ == "__main__":
    main()
