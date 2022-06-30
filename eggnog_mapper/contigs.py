"""
Tools for handling contigs and contig files
"""
from Bio import SeqIO

def contig_length_generator(input_fasta, min_contig_length):
    """
    :param input_fasta: FASTA file path
    :param min_contig_length: int minimum contig length to keep
    :returns:
      contig_filter_list - list of contig ids that pass the filter
      summary_file - path to contig length summary file
      contig_lengths - dict of contig id -> length (int), including filtered contigs
    """
    contig_filter_list = []
    fasta_file = open(input_fasta, 'r')
    contig_lengths = {}
    summary_file = f"{input_fasta.split('/')[-1].split('_')[0]}_contig_length_summary.tsv"
    with open(os.getcwd() + '/' + summary_file, 'w') as f:
        for rec in SeqIO.parse(fasta_file, 'fasta'):
            name = rec.id
            seq_len = len(rec)
            contig_lengths[name] = seq_len
            if seq_len < int(min_contig_length):
                contig_filter_list.append(name)
            f.write("{}\t{}\n".format(name, seq_len))
    fasta_file.close()
    return contig_filter_list, summary_file, contig_lengths


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
    Returns a list of contig ids that fall below the minimum contig length.
    """
    filtered_contigs = list()
    for contig_id, contig_length in contig_lengths.items():
        if contig_length < min_contig_length:
            filtered_contigs.append(contig_id)
    return filtered_contigs

