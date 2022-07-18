#!/usr/bin/env python

__author__ = "Sean P. Jungbluth"
__copyright__ = "Copyright 2022"
__maintainer__ = "Sean P. Jungbluth"
__email__ = "sjungbluth@lbl.com"

import os
import sys
import pandas
import argparse
import time
from collections import Counter, defaultdict
from .go_xref import GO_Xrefs
from .util import (
    isfloat,
    dict_to_csv
)
from .contigs import (
    load_contig_coverage,
    load_contig_filter,
    summarize_contig_lengths,
    get_allowed_contigs
)
from typing import List


### notes about eggNOG annotation output
# eggnog output categories
# query_name
# seed_eggNOG_ortholog
# seed_ortholog_evalue
# seed_ortholog_score
# best_tax_level
# Preferred_name
# GOs keep some connected to EC tails(!)
# EC keep - some connected to GO tails(!) painful
# KEGG_ko keep !
# KEGG_Pathway
# KEGG_Module keep
# KEGG_Reaction keep
# KEGG_rclass keep
# BRITE keep
# KEGG_TC keep
# CAZy keep
# BiGG_Reaction keep
# taxonomic_scope
# eggNOG_OGs keep
# best_eggNOG_OG
# COG_Functional_cat keep, but need to split more
# eggNOG_free_text_desc

## visually, sometimes seems like no tabs, but they ARE there...
#. e.g. try: grep "3.1.1.29" 1145081__M_70_70_06.annotations

### tool features:
# input annotation table
# weighted and unweighted coverage options
# different eggnog annotation cateogries

### todo
# add "skip contig taxa assignment" for speedup
# add tests e.g. https://realpython.com/python-testing/
# pca exploration of output
# add FAMA and other modules

def load_annotation_file(input_annotation_file: str) -> pandas.DataFrame:
    annotation_file_cols = [
        "query_name",
        "seed_eggNOG_ortholog",
        "seed_ortholog_evalue",
        "seed_ortholog_score",
        "best_tax_level",
        "Preferred_name",
        "GO",
        "EC",
        "KEGG_ko",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "KEGG_rclass",
        "BRITE",
        "KEGG_TC",
        "CAZy",
        "BiGG_Reaction",
        "taxonomic_scope",
        "eggNOG_OGs",
        "best_eggNOG_OG",
        "COG",
        "eggNOG_free_text_desc"
    ]
    return pandas.read_csv(os.path.abspath(input_annotation_file), delimiter='\t', header=None, names=annotation_file_cols)

# needed to detect transition to a new contig so that summary step executes to combine with coverage data on a per-contig basis
def _contig_summary_trigger(annotation_data, i):
    summarize_contig = False
    current_query_name = "_".join(annotation_data["query_name"][i].split("_")[:-1])  # store current query ID
    if i == len(annotation_data["query_name"]) - 1:  # detect if next ID end of list
        summarize_contig = True
    else:
        next_query_name = "_".join(annotation_data["query_name"][i+1].split("_")[:-1])
        if current_query_name != next_query_name:
            summarize_contig = True
    return(summarize_contig, current_query_name)

def _combine_and_add_dictionaries(dict_a, dict_b):
    dictcountA = Counter(dict_a)
    dictcountB = Counter(dict_b)
    dict_combined = dictcountA + dictcountB
    return dict_combined

def _determine_consensus_taxonomy(taxonomic_scope, taxonomy_consensus_threshold, contig_running_count):
    taxonomic_scope_output = dict(Counter(taxonomic_scope))
    consensus_taxonomy = "NA"
    for x in range(0, len(list(taxonomic_scope_output.keys()))):
        hit_frequency = list(taxonomic_scope_output.values())[x] / contig_running_count
        if hit_frequency >= taxonomy_consensus_threshold:
            consensus_taxonomy = list(taxonomic_scope_output.keys())[x]
            break
        consensus_taxonomy = "NA"
        hit_frequency = "NA"
    if isfloat(hit_frequency): # rarely, but can happen that there is no consensus taxonomy, so need to confirm
        hit_frequency = round(hit_frequency,2)
    return consensus_taxonomy, hit_frequency

# count features and multiple by coverage value, add to final summary dict
def _summarize_contig_module(f, summary_table, running_contig_non_hits, contig_running_count, cov_method, coverages, contig_id, summary_table_final_count, contig_length_data, taxonomic_scope, taxonomy_consensus_threshold):
    """
    :param f: output file pointer
    :param summary_table:
    running_contig_non_hits
    contig_running_count
    cov_method
    :param coverages: dict keys = contig ids, values = contig coverages
    contig_id
    summary_table_final_count
    :param contig_length_data: dict keys = contig ids, values = contig lengths
    taxonomic_scope
    taxonomy_consensus_threshold
    """
    summary_table_output = dict(Counter(summary_table))  # covert list to dictionary with feature counts
    feature_hit_freq = 1 - (running_contig_non_hits / contig_running_count)  # calculate frequency of genes with features on a contig
    contig_length = "NA"
    # determine the majority taxonmony for a contig
    consensus_taxonomy, consensus_taxonomy_frequency = _determine_consensus_taxonomy(taxonomic_scope, taxonomy_consensus_threshold, contig_running_count)
    if contig_id in contig_length_data:
        contig_length = contig_length_data[contig_id]
    if cov_method == "weighted":
        # multiply dictionary values by coverage for that contig
        if contig_id not in coverages:
            pass
        else:
            contig_coverage = coverages[contig_id]
            summary_table_output.update((x, y*contig_coverage) for x, y in summary_table_output.items())  # multiply dictionary values by contig coverage
            # need to add summary_table_output dict to new summary_table_final_count dict
            summary_table_final_count = _combine_and_add_dictionaries(summary_table_output, summary_table_final_count)
            f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(contig_id, contig_length, contig_coverage, round(feature_hit_freq,2), consensus_taxonomy, consensus_taxonomy_frequency))
    else:  # unweighted mode, no need to incorporate contig coverage information
        summary_table_final_count = _combine_and_add_dictionaries(summary_table, summary_table_final_count)
        f.write("{}\t{}\t{}\t{}\t{}\n".format(contig_id, contig_length, round(feature_hit_freq,2), consensus_taxonomy, consensus_taxonomy_frequency))
    # clear list of features for next contig
    summary_table = []
    taxonomic_scope = []
    contig_running_count = 0
    running_contig_non_hits = 0
    return summary_table_final_count, summary_table, contig_running_count, running_contig_non_hits, taxonomic_scope


# convert large list into dictionary and translate to pandas df for export
# need to generalize more - currently over-reliant on Torben metadone naming scheme
def _munge_prepare_and_export_count_tables(summary_table_final_count, input_annotation, input_cat, cov_method, binary_output, min_contig_length, min_contig_coverage):
    summary_table_final_count_dict = dict(Counter(summary_table_final_count))
    if binary_output == "Yes":
        for x in summary_table_final_count_dict:
            summary_table_final_count_dict[x] = 1
            dict_to_csv(summary_table_final_count_dict, f"{os.getcwd()}/{input_annotation.split('/')[-1].split('_')[0]}_len{min_contig_length}_cov{min_contig_coverage}_{input_cat}_binary_direct_count_table.csv")
    else:
        dict_to_csv(summary_table_final_count_dict, f"{os.getcwd()}/{input_annotation.split('/')[-1].split('_')[0]}_len{min_contig_length}_cov{min_contig_coverage}_{input_cat}_{cov_method}_direct_count_table.csv")

def _convert_go_count_table_to_other_annotation(summary_table_final_count, input_annotation, input_cat, cov_method, go_xref_tables, binary_output, min_contig_length, min_contig_coverage):
    print("\nGO count table production complete - converting GO to other name spaces using translation tables.")
    use_binary = binary_output == "Yes"
    for go_xref in go_xref_tables.keys():
        print(f"Converting...{go_xref}...")
        xref_summary_table_final_count_dict = defaultdict(int)
        xref_table = go_xref_tables[go_xref]

        for go_term in summary_table_final_count.keys():
            x_refs = xref_table[go_term]
            for x_ref in x_refs:
                xref_summary_table_final_count_dict[x_ref] += summary_table_final_count[go_term]

        outfile_path = f"{os.getcwd()}/{input_annotation.split('/')[-1].split('_')[0]}_len{min_contig_length}_cov{min_contig_coverage}_{input_cat}_{'binary' if use_binary else cov_method}_{go_xref}_count_table.csv"
        if use_binary:
            for x in xref_summary_table_final_count_dict:
                xref_summary_table_final_count_dict[x] = 1
        dict_to_csv(xref_summary_table_final_count_dict, outfile_path)

def scan_and_summarize_output(
    output_file_prefix,  # like "<mg_id>_len2000_cov5_GO_weighted"
    input_annotation,    # input annotation file name
    input_cat,           # input category - string
    cov_method,          # one of "weighted", "unweighted"
    annotation_data,     #
    coverages,
    contigs_allowed,     # set of used contigs
    make_go_xref,
    binary_output,
    contig_length_data,
    taxonomy_consensus_threshold,
    min_contig_length,
    min_contig_coverage,
    go_xrefs
):

    # prep starting variables
    summary_table = []
    taxonomic_scope = []
    summary_table_final_count = {}
    contig_running_count = 0
    running_contig_non_hits = 0

    cols = ["contig_id", "contig_length", "feature_hit_freq", "consensus_taxonomy", "consensus_taxonomy_frequency"]
    if cov_method == "weighted":
        cols.insert(2, "avg_contig_coverage")

    summary_file = f"{output_file_prefix}_summary.tsv"
    with open(summary_file, 'w') as f:
        f.write("\t".join(cols) + "\n")
        total_annotation_queries = len(annotation_data["query_name"])
        for i in range(total_annotation_queries):  # for each row in table
            if i % 1000 == 0:
                percent_done = round(((i+1)*100)/(total_annotation_queries),1)
                print("-----> Processing record {} of {}.....{} complete".format(i+1, total_annotation_queries, f"{percent_done} %"))
            summarize_contig, contig_id = _contig_summary_trigger(annotation_data, i)  # determine if last annotation of the current contig
            contig_running_count += 1
            if pandas.isna(annotation_data[input_cat][i]):  # if na, skip
                running_contig_non_hits += 1
            else:  # if not na, skip
                annotations = annotation_data[input_cat][i].split(",")
                for anno in annotations:
                    if input_cat == "COG":
                        if len(anno) > 1:
                            summary_table.append(anno)
                    else: # GO and other categories
                        summary_table.append(anno)
            taxonomic_scope.append(annotation_data["taxonomic_scope"][i])
            if summarize_contig:  # if end of contig reached, count features and multiple by coverage value, add to final summary dict
                if contig_id in contigs_allowed:
                    (summary_table_final_count,
                    summary_table,
                    contig_running_count,
                    running_contig_non_hits,
                    taxonomic_scope) = _summarize_contig_module(f, summary_table, running_contig_non_hits, contig_running_count, cov_method, coverages, contig_id, summary_table_final_count, contig_length_data, taxonomic_scope, taxonomy_consensus_threshold)
                else: # skip contig and reset relevant variables
                    taxonomic_scope = []
                    summary_table = []
                    contig_running_count = 0
                    running_contig_non_hits = 0
        _munge_prepare_and_export_count_tables(summary_table_final_count, input_annotation, input_cat, cov_method, binary_output, min_contig_length, min_contig_coverage)
        if make_go_xref == "Yes" and input_cat == "GO":
            _convert_go_count_table_to_other_annotation(summary_table_final_count, input_annotation, input_cat, cov_method, go_xrefs, binary_output, min_contig_length, min_contig_coverage)

def run_mapper(args):
    starttime = time.time()
    if args.binary_output == "Yes" and args.use_coverage:
        print("Binary output selected, adjusted coverage setting to maximize program performance")
        args.use_coverage = False
    if (not args.input_cat == "ALL_CATEGORIES") and (not args.input_cat == "GO") and (not args.input_cat == "EC") and (not args.input_cat == "KEGG_ko") and (not args.input_cat == "COG") and (not args.input_cat == "KEGG_Module") and (not args.input_cat == "KEGG_Reaction") and (not args.input_cat == "KEGG_rclass") and (not args.input_cat == "BRITE") and (not args.input_cat == "KEGG_TC") and (not args.input_cat == "CAZy") and (not args.input_cat == "BiGG_Reaction") and (not args.input_cat == "eggNOG_OGs"):
        sys.exit("Bad category selection, please review options. Exiting...")
    input_fasta = args.input_contig_fasta
    print('''
***************************************************************************
*                                App start                                *
***************************************************************************
    ''')
    print("Parameters used to run App:")
    print("")
    print("    Input annotation: " + str(args.input_annotation))
    print("    Input coverage: " + str(args.input_coverage))
    print("    Input contig_fasta: " + str(input_fasta))
    print("    Input min contig length: " + str(args.min_contig_length))
    print("    Input min contig coverage: " + str(args.min_contig_coverage))
    print("    Input contig filter: " + str(args.contig_filter))
    print("    Input eggnog mapper category: " + str(args.input_cat))
    print("    Use contig coverage: " + str(args.use_coverage))
    print("    Binary output: " + str(args.binary_output))
    print("    Contig taxonomy consensus threshold: " + str(args.taxonomy_consensus_threshold))
    print("    Make GO cross-ref tables: " + str(args.make_go_xref))
    print("")

    output_dir = os.getcwd()

    # contig length calculation
    input_fasta_prefix = os.path.basename(input_fasta).split("_")[0]
    fasta_summary_file = os.path.join(f"{output_dir}", f"{input_fasta_prefix}_contig_length_summary.tsv")

    contig_lengths = summarize_contig_lengths(input_fasta, fasta_summary_file)
    annotation_data = load_annotation_file(args.input_annotation)
    coverage_data = dict()
    if args.use_coverage or args.min_contig_coverage > 0:
        coverage_data = load_contig_coverage(args.input_coverage)
    filtered_contigs = set()
    if args.contig_filter is not None:
        filtered_contigs = load_contig_filter(args.contig_filter)
    contigs_allowed = get_allowed_contigs(contig_lengths, coverage_data, filtered_contigs, args.min_contig_coverage, args.min_contig_length)

    input_category_options = [
        "GO", "EC", "KEGG_ko", "COG", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "eggNOG_OGs"
    ]

    go_xrefs = GO_Xrefs()
    cov_method = "weighted" if args.use_coverage else "unweighted"

    # run main function to summarize eggnog mapper annotation output on a per contig basis
    if args.input_cat == "ALL_CATEGORIES":
        cat_options = input_category_options
    else:
        cat_options = [args.input_cat]
    for cat in cat_options:
        args.input_cat = cat
        print(f"Generating {args.input_cat} table(s) from eggNOG-mapper data")
        # import Gene Ontology cross-reference tables
        if args.input_cat == "GO" and args.make_go_xref != "No":
            go_xref_table_list = go_xrefs.xrefs
        else:
            go_xref_table_list = ""

        output_file_prefix = f"{os.getcwd()}/{args.input_annotation.split('/')[-1].split('_')[0]}_len{args.min_contig_length}_cov{args.min_contig_coverage}_{args.input_cat}_{cov_method}"

        scan_and_summarize_output(
            output_file_prefix,
            args.input_annotation,
            args.input_cat,
            cov_method,
            annotation_data,
            coverage_data,
            contigs_allowed,
            args.make_go_xref,
            args.binary_output,
            contig_lengths,
            args.taxonomy_consensus_threshold,
            args.min_contig_length,
            args.min_contig_coverage,
            go_xref_table_list
        )
    print('''
***************************************************************************
*                                 App end                                 *
***************************************************************************
    ''')
    print("App ran in %s seconds." % round((time.time() - starttime), 2))
    return 0

def get_mapper_args(arg_list: List[str]):
    if arg_list is None:
        arg_list = []
    description_text = "eggnog_to_feature_table converts eggnog mapper output to a feature count table, with or without weighting using contig abundance values."
    description_text += "\n" + "-"*len(description_text)
    parser = argparse.ArgumentParser(
        prog="eggnog_to_feature_table",
        usage=("%(prog)s.py \ \n"
            "\t-i [input_annotation] \ \n"
            "\t--input_cov [input_coverage] \ \n"
            "\t--min_contig_length [min_contig_length] \ \n"
            "\t--min_contig_coverage [min_contig_coverage] \ \n"
            "\t--contig_filter [contig_filter] \ \n"
            "\t--eggnog_catgory [input_cat] \ \n"
            "\t--use_cov [use_coverage] \ \n"
            "\t--binary [binary_output] \ \n"
            "\t--contig_taxa_threshold [taxonomy_consensus_threshold] \ \n"
            "\t--go_xref [make_go_xref] \ \n"),
        description=description_text,
        formatter_class=argparse.RawTextHelpFormatter
    )

    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-i", dest="input_annotation", help="Indicate the eggnog mapping annotation file for input. ", required=False)
    parser.add_argument("--input_contig_fasta", dest="input_contig_fasta", help="Indicate the contig fasta file to calculate sequence length. ")
    parser.add_argument("--input_cov", dest="input_coverage", help="Indicate the tab-separted contig coverage file for input. ")
    parser.add_argument("--min_contig_length", dest="min_contig_length", type=int, default=2000, help="Indicate the minimum contig length to include in final count table (default: 2000).")
    parser.add_argument("--min_contig_coverage", dest="min_contig_coverage", type=int, default=5, help="Indicate the minimum average contig coverage required to include in final count table (default: 5).")
    parser.add_argument("--contig_filter", dest="contig_filter", default=None, help="Indicate the contigs to retain for the output (i.e. contigs not listed will be filtered before count tables produced).")
    parser.add_argument("--eggnog_category", dest="input_cat", default="GO", help="Indicate the eggnog mapper annotation ontology to use for mapping. (options: ALL_CATEGORIES, GO, EC, KEGG_ko, COG, KEGG_Module, KEGG_Reaction, KEGG_rclass, BRITE, KEGG_TC, CAZy, BiGG_Reaction, eggNOG_OGs (default: GO))")
    parser.add_argument("--use_cov", dest="use_coverage", default=False, action="store_true", help="Indicate if contig coverage information should be used to produce a weighted output count table. (default: False)")
    parser.add_argument("--binary", dest="binary_output", default="No", help="Indicate if the output should just indicate hit presence/absence. (default: No)")
    parser.add_argument("--contig_taxa_threshold", dest="taxonomy_consensus_threshold", default=0.5, help="Indicate the frequency threshold for determining contig consensus taxonomy. (default: 0.5)")
    parser.add_argument("--go_xref", dest="make_go_xref", default="Yes", help="Indicate if input_cat is GO then generate cross-reference tables. (default: Yes)")
    parser.add_argument("--go_xref_loc", dest="go_xref_loc", default=".", help="Path to a directory containing GO cross-reference tables")
    parser.add_argument("--version", action="version", version='%(prog)s v2.0')
    parser.add_argument("--verbose", dest="verbose", default=False, action="store_true", help="Extra verbose output")

    if len(arg_list) < 2:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args(arg_list)
    # make args a little more sensible
    if args.binary_output == "Yes":
        args.use_coverage = False
        print("Binary output is selected, not using coverage for weighting")
    return args

if __name__ == '__main__':
    sys.exit(run_mapper(get_mapper_args(sys.argv[1:])))
