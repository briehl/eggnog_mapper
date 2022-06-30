from typing import List

def unique_list(input: List[str]) -> List[str]:
    """
    Make a unique set of inputs.
    """
    return list(set(input))

def isfloat(num):
    """
    A little wrapper around float() that just returns
    True or False instead of throwing exceptions
    """
    try:
        float(num)
        return True
    except ValueError:
        return False

