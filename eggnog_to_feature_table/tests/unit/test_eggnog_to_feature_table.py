from copy import copy
import pytest
import os
from eggnog_to_feature_table.eggnog_to_feature_table import (
    get_mapper_args,
    run_mapper
)

DEFAULT_ARGS = {
    "input_annotation": None,
    "input_contig_fasta": None,
    "input_coverage": None,
    "min_contig_length": 2000,
    "min_contig_coverage": 5,
    "contig_filter": None,
    "input_cat": "GO",
    "use_coverage": False,
    "binary_output": "No",
    "taxonomy_consensus_threshold": 0.5,
    "make_go_xref": "Yes",
    "verbose": False
}

ARGS_CASES = [
    (["-i", "foo.tsv"], {"input_annotation": "foo.tsv"}),
    (["--binary", "Yes", "--use_cov"], {
        "use_coverage": False,
        "binary_output": "Yes"
    }),
    (["-i", "foo.tsv", "--use_cov"], {"input_annotation": "foo.tsv", "use_coverage": True})
]

cur_dir = os.path.dirname(__file__)

TEST_INPUT_ANNOTATION = os.path.join(cur_dir, "../data/test_annotation.tsv")
TEST_INPUT_FASTA = os.path.join(cur_dir, "../data/test_fasta.fasta")
TEST_INPUT_COVERAGE =os.path.join(cur_dir,  "../data/test_coverage.tsv")

@pytest.mark.parametrize("input,expected", ARGS_CASES)
def test_get_mapper_args(input, expected):
    expected_args = copy(DEFAULT_ARGS)
    for k, v in expected.items():
        expected_args[k] = v
    args = vars(get_mapper_args(input))
    for k, v in expected_args.items():
        assert args[k] == v

def test_run_mapper_weighted_coverage():
    inputs = ["-i", TEST_INPUT_ANNOTATION, "--input_contig_fasta", TEST_INPUT_FASTA, "--input_cov", TEST_INPUT_COVERAGE, "--use_cov", "--go_xref", "Yes", "--eggnog_category", "GO"]
    args = get_mapper_args(inputs)
    run_mapper(args)

def test_run_mapper_unweighted():
    pass
