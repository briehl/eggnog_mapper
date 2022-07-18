import pandas as pd
import os

def load_annotations(input_annotation_file: str) -> pd.DataFrame:
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
    return pd.read_csv(os.path.abspath(input_annotation_file), delimiter='\t', header=None, names=annotation_file_cols)

def query_to_contig_id(query: str) -> str:
    """
    Assumes a query_id with format 1059893_contig_1_1 where the real id is 1059893_contig_1.
    Effectively, strip the last instance of _<value>

    Just a convenience method, might be adaptable to other formats later.
    """
    split_query = query.split("_")
    return "_".join(split_query[:-1])
