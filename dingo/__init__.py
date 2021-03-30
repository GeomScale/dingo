from dingo.fva import slow_fva
from dingo.fba import slow_fba
from dingo.loading_models import read_json_file
from dingo.inner_ball import slow_inner_ball
from dingo.nullspace import nullspace_dense, nullspace_sparse
from dingo.scaling import gmscale, apply_scaling, remove_almost_redundant_facets
from dingo.parse import dingo_args

from volestipy import HPolytope


def main_pipeline():
    args = dingo_args()

    return args


if __name__ == "__main__":

    b = main_pipeline()
