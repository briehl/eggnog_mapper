from eggnog_to_feature_table.contigs import (
    summarize_contig_lengths,
    filter_contig_lengths
)
from os import path
import tempfile

TEST_FASTA_FILE = path.join(path.dirname(__file__), "..", "data", "simple_file.fasta")
TEST_CONTIG_LENS = {
    "foo": 20,
    "bar": 40,
    "baz": 60
}
TEST_OUTPUT_FILE = "_test_simple_contigs.tsv"

def test_summarize_contig_length():
    with tempfile.TemporaryDirectory() as tmpdirname:
        test_outfile = path.join(tmpdirname, TEST_OUTPUT_FILE)
        contig_lens = summarize_contig_lengths(TEST_FASTA_FILE, test_outfile)
        assert contig_lens == TEST_CONTIG_LENS
        with open(test_outfile, "r") as testfile:
            lines = testfile.readlines()
            assert len(lines) == 3
            test_file_lens = dict()
            for line in lines:
                [key, val] = line.split('\t')
                test_file_lens[key] = int(val)
            assert test_file_lens == TEST_CONTIG_LENS

def test_filter_contig_lengths():
    filtered = filter_contig_lengths(TEST_CONTIG_LENS, 60)
    assert filtered == set(["foo", "bar"])
