from typing import Any, Dict, Union
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

def dict_to_csv(d: Dict[str, Union[int, str, float]], outfile_path: str) -> None:
    with open(outfile_path, 'w') as outfile:
        for k, v in d.items():
            outfile.write(f"{k},{v}\n")
