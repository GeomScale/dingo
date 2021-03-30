import numpy as np
import sys
import os
from dingo.fva import slow_fva
from dingo.fba import slow_fba
from dingo.loading_models import read_json_file
from dingo.inner_ball import slow_inner_ball
from dingo.nullspace import nullspace_dense, nullspace_sparse
from dingo.scaling import (
    gmscale,
    apply_scaling,
    remove_almost_redundant_facets,
    map_samples_to_steady_states,
)
from dingo.parser import dingo_args

try:
    import gurobipy
    from dingo.gurobi_based_implementations import fast_fba, fast_fva, fast_inner_ball
except ImportError as e:
    pass

from volestipy import HPolytope


def main_pipeline(args):

    if args.solver == "gurobi":
        try:
            import gurobipy
        except ImportError:
            print("Library gurobi is not available.")
            sys.exit(1)

    if args.solver != "gurobi" and args.solver != "scipy":
        raise Exception("An unknown solver requested.")

    if args.nullspace != "sparseQR" and args.nullspace != "scipy":
        raise Exception("An unknown method to compute the nullspace requested.")

    metabolic_network = read_json_file(args.metabolic_network)

    try:
        lb = metabolic_network[0]
        ub = metabolic_network[1]
        S = metabolic_network[2]
        biomass_index = metabolic_network[5]
        biomass_function = metabolic_network[6]
    except LookupError:
        print("An unexpected error occured when reading the input file.")
        sys.exit(1)

    if args.solver == "scipy":
        fva_res = slow_fva(lb, ub, S, biomass_function)
    elif args.solver == "gurobi":
        fva_res = fast_fva(lb, ub, S, biomass_function)

    A = fva_res[0]
    b = fva_res[1]
    Aeq = fva_res[2]
    beq = fva_res[3]
    min_fluxes = fva_res[4]
    max_fluxes = fva_res[5]

    if A.shape[0] != b.size or A.shape[1] != Aeq.shape[1] or Aeq.shape[0] != beq.size:
        raise Exception("FVA failed.")

    if args.nullspace == "sparseQR":
        nullspace_res = nullspace_sparse(Aeq, beq)
    elif args.nullspace == "scipy":
        nullspace_res = nullspace_dense(Aeq, beq)

    N = nullspace_res[0]
    N_shift = nullspace_res[1]

    if A.shape[1] != N.shape[0] or N.shape[0] != N_shift.size or N.shape[1] <= 1:
        raise Exception(
            "The computation of the matrix of the right nullspace of the stoichiometric matrix failed."
        )

    product = np.dot(A, N_shift)
    b = np.subtract(b, product)
    A = np.dot(A, N)

    res = remove_almost_redundant_facets(A, b)
    A = res[0]
    b = res[1]

    try:
        res = gmscale(A, 0.99)
        res = apply_scaling(A, b, res[0], res[1])
        A = res[0]
        b = res[1]
        N = np.dot(N, res[2])

        res = remove_almost_redundant_facets(A, b)
        A = res[0]
        b = res[1]
    except:
        print("gmscale failed to compute a good scaling.")

    p = HPolytope(A, b)

    if args.solver == "scipy":
        A, b, T, T_shift, samples = p.slow_mmcs(
            args.effective_sample_size, args.psrf_check
        )
    elif args.solver == "gurobi":
        A, b, T, T_shift, samples = p.fast_mmcs(
            args.effective_sample_size, args.psrf_check
        )

    if (
        A.shape[1] != N.shape[1]
        or A.shape[0] != b.size
        or A.shape[1] != T.shape[0]
        or T.shape[1] != T_shift.size
        or samples.shape[0] != A.shape[1]
    ):
        raise Exception(
            "An unexpected error occured in Multiphase Monte Carlo Sampling algorithm."
        )

    steady_states = map_samples_to_steady_states(samples, T, T_shift, N, N_shift)

    return A, b, T, T_shift, N, N_shift, min_fluxes, max_fluxes, samples, steady_states


def dingo_main():

    args = dingo_args()
    # print(args)

    result_obj = main_pipeline(args)

    if args.output_directory == None:
        output_path_dir = os.getcwd()
    else:
        output_path_dir = args.output_directory

    if os.path.isdir(output_path_dir) == False:
        os.mkdir(output_path_dir)

    # Move to the model's output directory
    os.chdir(output_path_dir)

    # np.save(os.path.join(output_path_dir,'dingo_output.npy'), result_obj)


if __name__ == "__main__":

    dingo_main()
