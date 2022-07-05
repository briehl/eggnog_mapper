from typing import Any
from math import isnan

def isfloat(num: Any) -> bool:
    """
    A little wrapper around float() that just returns
    True or False instead of throwing exceptions
    """
    try:
        float_num = float(num)
        return not isnan(float_num)
    except (ValueError, TypeError):
        return False
