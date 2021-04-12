# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import sys
from dingo.fva import slow_fva
from dingo.loading_models import read_json_file, read_mat_file
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


def from_model_to_steady_states_pipeline(args):
    """A function to compute the full dimensional polytope from a model (metabolic network), sample
    from it and to transform the sample to steady states using (M)ultiphase (M)onte (C)arlo (S)ampling algorithm

    Keyword arguments:
    args -- a Namespace that contains the input file and the requested methods, solvers and MMCS parameters; see parser.py

    Keyword outputs:
    A -- the matrix of the full dimensional polytope, s.t. Ax <= b
    b -- the vector of the full dimensional polytope, s.t. Ax <= b
    N -- the matrix of the right nullspace of the stoichiometric matrix S
    N_shift -- a solution of the linear system defined by the augmented stoichiometric matrix
    T -- the matrix of the linear transformation that maps the full dimensional polytope to a rounded polytope
    T_shift -- the shifting vector of the linear transformation that maps the full dimensional polytope to a rounded polytope
    min_fluxes -- the minimum values of the fluxes of each reaction
    max_fluxes -- the maximum values of the fluxes of each reaction
    steady_states -- a matrix that contains column-wise the generated steady states
    """

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

    if args.metabolic_network[-3:] != "mat" and args.metabolic_network[-4:] != "json":
        raise Exception("An unknown format file given.")

    if args.metabolic_network[-4:] == "json":
        metabolic_network = read_json_file(args.metabolic_network)
    else:
        metabolic_network = read_mat_file(args.metabolic_network)

    try:
        lb = metabolic_network[0]
        ub = metabolic_network[1]
        S = metabolic_network[2]
        metabolites = metabolic_network[3]
        reactions = metabolic_network[4]
        biomass_index = metabolic_network[5]
        biomass_function = metabolic_network[6]
    except LookupError:
        print("An unexpected error occured when reading the input file.")
        sys.exit(1)

    if args.unbiased_analysis:
        biomass_function = np.zeros(S.shape[1])

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

    if args.preprocess_only:
        return A, b, N, N_shift, min_fluxes, max_fluxes, metabolites, reactions

    p = HPolytope(A, b)

    if args.solver == "scipy":
        A, b, T, T_shift, samples = p.slow_mmcs(
            int(args.effective_sample_size),
            args.psrf_check,
            args.parallel_mmcs,
            int(args.num_threads),
        )
    elif args.solver == "gurobi":
        A, b, T, T_shift, samples = p.fast_mmcs(
            int(args.effective_sample_size),
            args.psrf_check,
            args.parallel_mmcs,
            int(args.num_threads),
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

    steady_states = map_samples_to_steady_states(samples, N, N_shift)

    return A, b, N, N_shift, T, T_shift, samples, min_fluxes, max_fluxes, metabolites, reactions, steady_states


def from_polytope_to_steady_states_pipeline(
    args, A, b, N, N_shift, T=None, T_shift=None
):
    """A function to sample from the fulldimensional polytope derived from a model (metabolic network)
    and to transform the sample to steady states using (M)ultiphase (M)onte (C)arlo (S)ampling algorithm

    Keyword arguments:
    args -- a Namespace that contains the input file and the requested methods, solvers and MMCS parameters; see parser.py
    A -- the matrix of the full dimensional polytope, s.t. Ax <= b
    b -- the vector of the full dimensional polytope, s.t. Ax <= b
    N -- the matrix of the right nullspace of the stoichiometric matrix S
    N_shift -- a solution of the linear system defined by the augmented stoichiometric matrix
    T -- optional; the matrix of the linear transformation that maps the full dimensional polytope to a rounded polytope
    T_shift -- optional; the shifting vector of the linear transformation that maps the full dimensional polytope to a rounded polytope

    Keyword outputs:
    A -- the matrix of the rounded full dimensional polytope, s.t. Ax <= b
    b -- the vector of the rounded full dimensional polytope, s.t. Ax <= b
    N -- the matrix of the right nullspace of the stoichiometric matrix S
    N_shift -- a solution of the linear system defined by the augmented stoichiometric matrix
    T -- the matrix of the linear transformation that maps the full dimensional polytope to the input rounded polytope
    T_shift -- the shifting vector of the linear transformation that maps the full dimensional polytope to the input rounded polytope
    steady_states -- a matrix that contains column-wise the generated steady states
    """

    if args.solver == "gurobi":
        try:
            import gurobipy
        except ImportError:
            print("Library gurobi is not available.")
            sys.exit(1)

    if args.solver != "gurobi" and args.solver != "scipy":
        raise Exception("An unknown solver requested.")

    p = HPolytope(A, b)

    if args.solver == "scipy":
        A, b, Tr, Tr_shift, samples = p.slow_mmcs(
            int(args.effective_sample_size),
            args.psrf_check,
            args.parallel_mmcs,
            int(args.num_threads),
        )
    elif args.solver == "gurobi":
        A, b, Tr, Tr_shift, samples = p.fast_mmcs(
            int(args.effective_sample_size),
            args.psrf_check,
            args.parallel_mmcs,
            int(args.num_threads),
        )

    if T is not None and T_shift is not None:
        Tr = np.dot(T, Tr)
        Tr_shift = T_shift + Tr_shift

    if (
        A.shape[1] != N.shape[1]
        or A.shape[0] != b.size
        or A.shape[1] != Tr.shape[0]
        or Tr.shape[1] != Tr_shift.size
        or samples.shape[0] != A.shape[1]
    ):
        raise Exception(
            "An unexpected error occured in Multiphase Monte Carlo Sampling algorithm."
        )

    steady_states = map_samples_to_steady_states(samples, N, N_shift, T, T_shift)

    return A, b, N, N_shift, Tr, Tr_shift, samples, steady_states


def fva_pipeline(args):

    if args.solver == "gurobi":
        try:
            import gurobipy
        except ImportError:
            print("Library gurobi is not available.")
            sys.exit(1)

    if args.solver != "gurobi" and args.solver != "scipy":
        raise Exception("An unknown solver requested.")

    if args.metabolic_network[-3:] != "mat" and args.metabolic_network[-4:] != "json":
        raise Exception("An unknown format file given.")

    if args.metabolic_network[-4:] == "json":
        metabolic_network = read_json_file(args.metabolic_network)
    else:
        metabolic_network = read_mat_file(args.metabolic_network)

    try:
        lb = metabolic_network[0]
        ub = metabolic_network[1]
        S = metabolic_network[2]
        biomass_function = metabolic_network[6]
    except LookupError:
        print("An unexpected error occured when reading the input file.")
        sys.exit(1)

    if args.solver == "scipy":
        fva_res = slow_fva(lb, ub, S, biomass_function, args.opt_percentage)
    elif args.solver == "gurobi":
        fva_res = fast_fva(lb, ub, S, biomass_function, args.opt_percentage)

    return fva_res


def fba_pipeline(args):

    if args.solver == "gurobi":
        try:
            import gurobipy
        except ImportError:
            print("Library gurobi is not available.")
            sys.exit(1)

    if args.solver != "gurobi" and args.solver != "scipy":
        raise Exception("An unknown solver requested.")

    if args.metabolic_network[-3:] != "mat" and args.metabolic_network[-4:] != "json":
        raise Exception("An unknown format file given.")

    if args.metabolic_network[-4:] == "json":
        metabolic_network = read_json_file(args.metabolic_network)
    else:
        metabolic_network = read_mat_file(args.metabolic_network)

    try:
        lb = metabolic_network[0]
        ub = metabolic_network[1]
        S = metabolic_network[2]
        biomass_function = metabolic_network[6]
    except LookupError:
        print("An unexpected error occured when reading the input file.")
        sys.exit(1)

    if args.solver == "scipy":
        fba_res = fast_fba(lb, ub, S, biomass_function)
    elif args.solver == "gurobi":
        fba_res = fast_fba(lb, ub, S, biomass_function)

    return fba_res
