#!/usr/bin/env python3
"""
Universal script to process any virus family:
1. Extract top 1 hit from existing TSV
2. Re-run BLAST to get alignment coordinates
3. Download annotated GenBank reference genomes
4. Group contigs by reference genome
5. Convert GenBank files to FASTA format
"""

import os
import sys
import time
import glob
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

# Set your email for NCBI Entrez (required by NCBI)
Entrez.email = "user@institute.edu"

def extract_top_hits(tsv_file, output_file, top_n=1):
    """Extract top N hits per contig"""
    df = pd.read_csv(tsv_file, sep='\t')
    top_hits = df.groupby('query_id').head(top_n).reset_index(drop=True)
    top_hits.to_csv(output_file, sep='\t', index=False)

    print(f"  Total hits: {len(df)}")
    print(f"  Unique contigs: {df['query_id'].nunique()}")
    print(f"  Top {top_n} hits extracted: {len(top_hits)}")
    return top_hits

def blast_with_coordinates(seq_record, max_hits=1):
    """Run BLASTN and extract alignment coordinates"""
    print(f"  Running BLASTN for {seq_record.id} (length: {len(seq_record.seq)} bp)...")

    try:
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            sequence=str(seq_record.seq),
            hitlist_size=max_hits,
            megablast=True
        )

        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)

        hits = []
        for alignment in blast_record.alignments[:max_hits]:
            hsp = alignment.hsps[0]

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
        print(f"  Error: {e}")
        return []

def download_genbank_record(accession, output_dir):
    """Download GenBank file"""
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{accession}.gb")

    if os.path.exists(output_file):
        print(f"  {accession} already downloaded")
        return output_file

    try:
        print(f"  Downloading {accession}...")
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        genbank_data = handle.read()
        handle.close()

        with open(output_file, 'w') as f:
            f.write(genbank_data)

        time.sleep(0.5)
        return output_file

    except Exception as e:
        print(f"  Error downloading {accession}: {e}")
        return None

def check_contig_overlap(contig1, contig2):
    """Check if two contigs overlap on reference"""
    c1_start = min(contig1['subject_start'], contig1['subject_end'])
    c1_end = max(contig1['subject_start'], contig1['subject_end'])
    c2_start = min(contig2['subject_start'], contig2['subject_end'])
    c2_end = max(contig2['subject_start'], contig2['subject_end'])

    overlap = not (c1_end < c2_start or c2_end < c1_start)

    if overlap:
        overlap_start = max(c1_start, c2_start)
        overlap_end = min(c1_end, c2_end)
        overlap_length = overlap_end - overlap_start
        return True, overlap_length
    else:
        return False, 0

def convert_genbank_to_fasta(genbank_file, output_dir):
    """Convert GenBank to FASTA"""
    try:
        record = SeqIO.read(genbank_file, "genbank")
        basename = os.path.basename(genbank_file).replace('.gb', '')
        output_file = os.path.join(output_dir, f"{basename}_genome.fasta")
        SeqIO.write(record, output_file, "fasta")
        return output_file
    except Exception as e:
        print(f"  Error converting {genbank_file}: {e}")
        return None

def process_virus_family(virus_family):
    """Main processing function for a virus family"""

    # Set up paths
    base_dir = virus_family
    fasta_files = glob.glob(os.path.join(base_dir, "*.fas")) + \
                  glob.glob(os.path.join(base_dir, "*.fas.txt"))
    tsv_files = glob.glob(os.path.join(base_dir, "*.tsv"))

    if not fasta_files or not tsv_files:
        print(f"Error: Could not find FASTA or TSV files in {base_dir}")
        return

    fasta_file = fasta_files[0]
    tsv_file = tsv_files[0]

    print("=" * 70)
    print(f"PROCESSING: {virus_family}")
    print("=" * 70)
    print(f"FASTA: {fasta_file}")
    print(f"TSV: {tsv_file}")
    print()

    # Step 1: Extract top 1 hit (only if TSV has correct format)
    print("STEP 1: Checking existing TSV format")
    print("-" * 70)
    try:
        df_check = pd.read_csv(tsv_file, sep='\t', nrows=1)
        if 'query_id' in df_check.columns:
            print("  TSV has query_id column - extracting top 1 hits")
            top1_file = os.path.join(base_dir, f"{virus_family}_top1_hits.tsv")
            top1_df = extract_top_hits(tsv_file, top1_file, top_n=1)
        else:
            print("  TSV has different format - will run BLAST directly")
    except Exception as e:
        print(f"  Could not read TSV: {e}")
        print("  Will run BLAST directly")
    print()

    # Step 2: Re-run BLAST with coordinates
    print("STEP 2: Re-running BLAST to get alignment coordinates")
    print("-" * 70)
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"  Found {len(sequences)} sequences")

    all_results = []
    for i, seq_record in enumerate(sequences, 1):
        print(f"  [{i}/{len(sequences)}] {seq_record.id}...")
        hits = blast_with_coordinates(seq_record, max_hits=1)
        all_results.extend(hits)

        if i < len(sequences):
            print(f"  Waiting 10 seconds...")
            time.sleep(10)

    coords_file = os.path.join(base_dir, f"{virus_family}_with_coordinates.tsv")
    df = pd.DataFrame(all_results)
    df.to_csv(coords_file, sep='\t', index=False)
    print(f"  Saved {len(all_results)} results to {coords_file}")
    print()

    # Step 3: Download GenBank files
    print("STEP 3: Downloading annotated reference genomes")
    print("-" * 70)
    unique_refs = df['hit_accession'].unique()
    print(f"  Found {len(unique_refs)} unique reference genomes")

    ref_dir = os.path.join(base_dir, f"{virus_family}_reference_genomes")
    for accession in unique_refs:
        download_genbank_record(accession, ref_dir)
    print()

    # Step 4: Group contigs
    print("STEP 4: Grouping contigs by reference genome")
    print("-" * 70)
    grouped = df.groupby('hit_accession')

    for ref_acc, group in grouped:
        print(f"\n  Reference: {ref_acc}")
        print(f"    {group.iloc[0]['hit_description'][:80]}")
        print(f"    Number of contigs: {len(group)}")

        if len(group) > 1:
            print(f"    Contigs: {', '.join(group['query_id'].tolist())}")

            contigs_list = group.to_dict('records')
            print("    Checking overlaps:")
            for i, c1 in enumerate(contigs_list):
                for j, c2 in enumerate(contigs_list[i+1:], i+1):
                    overlap, overlap_len = check_contig_overlap(c1, c2)
                    if overlap:
                        print(f"      ⚠️  {c1['query_id']} and {c2['query_id']} OVERLAP by {overlap_len} bp")
                    else:
                        gap = min(abs(c1['subject_start'] - c2['subject_end']),
                                abs(c2['subject_start'] - c1['subject_end']))
                        print(f"      ✓  {c1['query_id']} and {c2['query_id']} separate (gap: {gap} bp)")

            print("    Alignment positions:")
            for _, row in group.iterrows():
                print(f"      {row['query_id']:15} -> ref:{row['subject_start']:6}-{row['subject_end']:6} ({row['strand']})")
        else:
            print(f"    Contig: {group.iloc[0]['query_id']}")
    print()

    # Step 5: Convert GenBank to FASTA
    print("STEP 5: Converting GenBank files to FASTA")
    print("-" * 70)
    genbank_files = glob.glob(os.path.join(ref_dir, "*.gb"))
    genome_dir = os.path.join(ref_dir, "genome_sequences")
    os.makedirs(genome_dir, exist_ok=True)

    for gb_file in genbank_files:
        accession = os.path.basename(gb_file).replace('.gb', '')
        fasta_out = convert_genbank_to_fasta(gb_file, genome_dir)
        if fasta_out:
            print(f"  ✓ {accession} -> {os.path.basename(fasta_out)}")

    print()
    print("=" * 70)
    print(f"COMPLETED: {virus_family}")
    print("=" * 70)
    print()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python process_virus_family.py <virus_family_name>")
        sys.exit(1)

    virus_family = sys.argv[1]
    process_virus_family(virus_family)
