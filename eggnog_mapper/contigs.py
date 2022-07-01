"""
Tools for handling contigs and contig files
"""
from Bio import SeqIO
import os

def summarize_contig_lengths(input_fasta_path, output_path):
    fasta_file = open(input_fasta_path, 'r')
    contig_lengths = {}

    with open(output_path, "w") as f:
        for rec in SeqIO.parse(fasta_file, 'fasta'):
            name = rec.id
            seq_len = len(rec)
            contig_lengths[name] = seq_len
            f.write("{}\t{}\n".format(name, seq_len))
    fasta_file.close()
    return contig_lengths

def filter_contig_lengths(contig_lengths, min_contig_length):
    """
    Returns a set of contig ids that fall below the minimum contig length.
    """
    filtered_contigs = set()
    for contig_id, contig_length in contig_lengths.items():
        if contig_length < min_contig_length:
            filtered_contigs.add(contig_id)
    return filtered_contigs

