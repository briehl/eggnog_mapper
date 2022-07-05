import re
from eggnog_mapper.go_xref import (
    GO_Xrefs,
    XREF_TYPES
)

def test_go_xrefs():
    # make sure it loads, has the relevant keys, and each maps from
    # a GO id string to a list of strings.
    xrefs = GO_Xrefs()
    assert len(xrefs.xrefs.keys()) == len(XREF_TYPES)
    for xref in XREF_TYPES:
        assert xrefs.xrefs[xref]
        for go_id, x_ref_list in xrefs.xrefs[xref].items():
            assert re.match("^GO:\\d+$", go_id)
            assert isinstance(x_ref_list, list)
            for term in x_ref_list:
                assert isinstance(term, str)
