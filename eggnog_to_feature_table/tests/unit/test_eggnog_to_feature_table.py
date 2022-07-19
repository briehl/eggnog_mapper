from contextlib import redirect_stdout
from copy import copy
import pytest
from os import path
from io import StringIO
from eggnog_to_feature_table.eggnog_to_feature_table import (
    get_mapper_args,
    run_mapper,
    build_data_summary
)
import tempfile
import json

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

cur_dir = path.dirname(__file__)

TEST_INPUT_ANNOTATION = path.join(cur_dir, "../data/test_annotation.tsv")
TEST_INPUT_FASTA = path.join(cur_dir, "../data/test_fasta.fasta")
TEST_INPUT_COVERAGE =path.join(cur_dir,  "../data/test_coverage.tsv")

SUMMARY_DATA_DIR = path.join(cur_dir, "../data/summary_cases")

SUMMARY_ANNOTATION = path.join(cur_dir, "../data/summary_annotation.tsv")
SUMMARY_FASTA = path.join(cur_dir, "../data/summary_fasta.fasta")
SUMMARY_COVERAGE =path.join(cur_dir,  "../data/summary_coverage.tsv")

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

def test_run_app_summary_only():
    """
    run the whole app with parameters, etc., expect only the summary mode to run
    """
    pass

@pytest.mark.parametrize("summary_case", ["simple", "missing_annotations", "missing_coverage", "missing_fasta", "no_overlap"])
def test_build_data_summary(summary_case):
    annotation_file = path.join(SUMMARY_DATA_DIR, summary_case, "annotation.tsv")
    fasta_file = path.join(SUMMARY_DATA_DIR, summary_case, "contigs.fasta")
    coverage_file = path.join(SUMMARY_DATA_DIR, summary_case, "coverage.tsv")
    expected_log_file = path.join(SUMMARY_DATA_DIR, summary_case, "expected_log.txt")
    expected_summary_file = path.join(SUMMARY_DATA_DIR, summary_case, "summary.tsv")
    expected_json_file = path.join(SUMMARY_DATA_DIR, summary_case, "result.json")

    with tempfile.TemporaryDirectory() as tmpdirname:
        summary_file_path = path.join(tmpdirname, "test_summary_file.tsv")
        summary_inputs = [
            "-i", annotation_file,
            "--input_contig_fasta", fasta_file,
            "--input_cov", coverage_file,
            "--min_contig_coverage", "5",
            "--min_contig_length", "10",
            "--summary_only",
            "--summary_file", summary_file_path
        ]
        with StringIO() as test_stdout, redirect_stdout(test_stdout):
            args = get_mapper_args(summary_inputs)
            summary = build_data_summary(args)
            summary_out = test_stdout.getvalue()

        # check the output file is as expected
        with open(summary_file_path, "r") as summary_file:
            summary_lines = summary_file.readlines()
        with open(expected_summary_file, "r") as expected_summary:
            expected_summary_lines = expected_summary.readlines()
        assert summary_lines == expected_summary_lines

        # check the result is as expected
        with open(expected_json_file) as expected_json:
            loaded_json = json.load(expected_json)
            assert loaded_json == summary

        summary_out_lines = summary_out.split("\n")
        print(summary_out_lines)
        print(len(summary_out_lines))

        assert summary_out_lines[0:3] == [
            f"Annotation file: {annotation_file}",
            f"FASTA file: {fasta_file}",
            f"Coverage file: {coverage_file}"
        ]
        # check the log file is as expected
        with open(expected_log_file) as expected_log:
            expected_log_lines = [s.strip() for s in expected_log.readlines()]
            assert expected_log_lines == summary_out_lines[3:]
