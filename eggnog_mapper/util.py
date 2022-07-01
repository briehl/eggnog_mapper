from typing import List

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

