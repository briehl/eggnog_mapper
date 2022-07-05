from eggnog_mapper.util import isfloat

def test_isfloat_false():
    for test_case in ["foo", "", "NaN", list(), dict(), None]:
        assert isfloat(test_case) == False

def test_isfloat_true():
    for test_case in [1.1, 1, 0, 0.00000001, -1.12345, "1", "1.1", "-1"]:
        assert isfloat(test_case) == True
