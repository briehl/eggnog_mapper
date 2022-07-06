import sys
from .eggnog_to_feature_table import (
    run_mapper,
    get_mapper_args
)

if __name__ == "__main__":
    args = get_mapper_args(sys.argv[1:])
    run_mapper(get_mapper_args(sys.argv[1:]))
