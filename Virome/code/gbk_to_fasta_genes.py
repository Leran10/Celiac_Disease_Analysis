from Bio import SeqIO
from Bio.Seq import Seq
import argparse

def get_position(pos):
    """Helper function to get position that works across BioPython versions"""
    try:
        return pos.position
    except AttributeError:
        return int(pos)

def extract_sequences(gbk_file, fasta_file, output_file):
    """
    Extract CDS sequences from GBK file and cross-reference with FASTA.
    
    Parameters:
    gbk_file (str): Path to input GenBank file
    fasta_file (str): Path to reference FASTA file
    output_file (str): Path to output FASTA file
    """
    # Read the reference sequence
    reference_record = next(SeqIO.parse(fasta_file, "fasta"))
    reference_seq = str(reference_record.seq)
    
    # Store extracted sequences
    extracted_sequences = []
    
    # Parse GenBank file
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                try:
                    # Get CDS function/product and ID
                    cds_function = feature.qualifiers.get("product", ["Unknown"])[0]
                    cds_id = feature.qualifiers.get("ID", ["Unknown"])[0]
                    
                    # Get sequence coordinates using version-safe method
                    start = get_position(feature.location.start)
                    end = get_position(feature.location.end)
                    strand = feature.location.strand
                    
                    # Extract sequence from reference
                    sequence = reference_seq[start:end]
                    
                    # Reverse complement if on negative strand
                    if strand == -1:
                        sequence = str(Seq(sequence).reverse_complement())
                    
                    # Create sequence record with modified format (removed "ID=")
                    description = f"{cds_id} CDS={cds_function} location={start}..{end} strand={'forward' if strand == 1 else 'reverse'}"
                    extracted_sequences.append(f">{description}\n{sequence}\n")
                    
                except Exception as e:
                    print(f"Error processing CDS {cds_function}: {str(e)}")
    
    # Write sequences to output file
    with open(output_file, 'w') as f:
        f.writelines(extracted_sequences)
    
    print(f"Extracted {len(extracted_sequences)} CDS sequences to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract CDS sequences from GBK and FASTA files")
    parser.add_argument("gbk_file", help="Input GenBank file")
    parser.add_argument("fasta_file", help="Reference FASTA file")
    parser.add_argument("output_file", help="Output FASTA file")
    
    args = parser.parse_args()
    extract_sequences(args.gbk_file, args.fasta_file, args.output_file)

if __name__ == "__main__":
    main()