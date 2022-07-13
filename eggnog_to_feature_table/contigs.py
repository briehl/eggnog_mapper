"""
Tools for handling contigs and contig files
"""
from Bio import SeqIO
from typing import Any, Dict, Set, Union
import os

def summarize_contig_lengths(input_fasta_path: str, output_path: str) -> Dict[str, int]:
    """
    Given a FASTA file path, this opens the file, reads all FASTA strings, and returns
    a dictionary with the ids mapped to the length of the FASTA feature.

    e.g. a file:
    > foo
    ATTGGC
    > bar
    AAGGCCTT

    would return the dict:
    {
        foo: 6,
        bar: 8
    }

    Note that this assumes unique ids in the FASTA file.

    This also prints out the results to the given output_path file, overwriting it if
    already present.

    :param input_fasta_path: Path to the input FASTA file
    :param output_path: Path to the output file
    :returns: dict (string -> int) as described above.
    """
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

def filter_contig_lengths(contig_lengths: Dict[str, int], min_contig_length: int) -> Set[str]:
    """
    Returns a Set of contig ids that fall below the minimum contig length.
    These would, conceivably, be removed from further analyses.

    :param contig_lengths: dict mapping from contig id -> int
    :param min_contig_length: the minimum contig length to pass the filter
    :returns: Set of contig id strings
    """
    filtered_contigs = set()
    for contig_id, contig_length in contig_lengths.items():
        if contig_length < min_contig_length:
            filtered_contigs.add(contig_id)
    return filtered_contigs

def load_contig_coverage(contig_cov_file: str) -> Dict[str, int]:
    coverages = {}
    with open(os.path.abspath(contig_cov_file), "r") as cov_file:
        for line in cov_file:
            row = line.strip().split('\t')
            coverages[row[0]] = int(row[1])
    return coverages

def load_contig_filter(contig_filter_file: str) -> Set[str]:
    contigs = set()
    with open(os.path.abspath(contig_filter_file), "r") as filter_file:
        for line in filter_file:
            contigs.add(line.strip())
    return contigs

def get_allowed_contigs(contig_lengths: Dict[str, int], contig_coverage: Dict[str, int], filtered_contigs: Set[str], min_coverage: int, min_length: int) -> Set[str]:
    if len(contig_lengths) != len(contig_coverage):
        print(
            "WARNING: number of contigs with known lengths is not equal to number of contigs with coverage data\n",
            "  Ignoring length or coverage requirement if those contig ids are not present."
        )
        print("Contigs with length and not coverage:")
        print("\n".join(contig_lengths.keys() - contig_coverage.keys()))
        print("Contigs with coverage and not length:")
        print("\n".join(contig_coverage.keys() - contig_lengths.keys()))
    allowed = set()
    all_contigs = contig_lengths.keys() | contig_coverage.keys()
    for contig in all_contigs:
        if contig_lengths.get(contig, min_length) >= min_length and \
            contig_coverage.get(contig, min_coverage) and \
            contig not in filtered_contigs:
            allowed.add(contig)
    return allowed

