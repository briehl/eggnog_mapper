import os
from collections import defaultdict

XREF_TYPES = [
    "ec2go",
    "hamap2go",
    "interpro2go",
    "kegg_reaction2go",
    "metacyc2go",
    "pfam2go",
    "pirsf2go",
    "prosite2go",
    "reactome2go",
    "rfam2go",
    "rhea2go",
    "smart2go",
    "um-bbd_enzymeid2go",
    "um-bbd_pathwayid2go",
    "um-bbd_reactionid2go",
    "uniprotkb_kw2go",
    "uniprotkb_sl2go",
    "unirule2go"
]

class GO_Xrefs:
    """
    Imports and manages the mappings from the Gene Ontology to various other namespaces.
    """
    def __init__(self, xref_dir="../GO_xref/"):
        self.xrefs = self._import_go_xref_tables(xref_dir)

    def _read_xref_table(self, path: str) -> dict:
        xref_table = defaultdict(list)
        with open(path) as f:
            lines = f.readlines()
            for line in lines:
                row = line.strip().split('\t')
                # 3 here - XRef_Annotation and GO_desc
                xref_table[row[2]].append(row[0])
        return xref_table

    # import all GO cross-reference tables
    def _import_go_xref_tables(self, go_directory="../GO_xref/"):
        """
        Imports all GO cross-reference tables. The table names are given in XREF_TYPES, and each
        cross-ref file is expected to be in go_directory/<file>.clean. Each file is expected to have
        the format:
        Namespace ID \t GO name \t GO id

        When loaded, the keys are all GO ids, and values are lists of 1 or more ids of the original
        namespace. These results are all available under xrefs[file]

        e.g.:

        all_xrefs = GO_Xrefs()
        all_xrefs.xrefs["ec2go"]["GO:0008465"] = ["EC:1.1.1.29"]
        """
        xrefs = {}

        for file in XREF_TYPES:
            filepath = os.path.join(os.path.dirname(__file__), go_directory, f"{file}.clean")
            xrefs[f"{file}_data"] = self._read_xref_table(filepath)
        return xrefs

