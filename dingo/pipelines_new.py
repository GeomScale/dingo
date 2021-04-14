# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import sys
import matplotlib.pyplot as plt
from dingo.fva import slow_fva
from dingo.loading_models import read_json_file, read_mat_file
from dingo.inner_ball import slow_inner_ball
from dingo.nullspace import nullspace_dense, nullspace_sparse
from dingo.scaling import gmscale
from dingo.utils import (
    apply_scaling,
    remove_almost_redundant_facets,
    map_samples_to_steady_states,
    get_matrices_of_low_dim_polytope,
    get_matrices_of_full_dim_polytope,
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

    model = metabolic_network(args.metabolic_network)

    sampler = network_sampler(model)

    if args.preprocess_only:
        sampler.get_polytope()

        return model, sampler

    steady_states = sampler.generate_steady_states(int(args.effective_sample_size), args.psrf_check, args.parallel_mmcs, int(args.num_threads))

    return model, sampler, steady_states


def from_polytope_to_steady_states_pipeline(args, sampler):
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

    steady_states = sampler.generate_steady_states(int(args.effective_sample_size), args.psrf_check, args.parallel_mmcs, int(args.num_threads))

    return steady_states

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

    model = metabolic_network(args.metabolic_network)

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


def plot_histogram(reaction_fluxes, reaction, n_bins):

    plt.figure(figsize=(7, 7))

    n, bins, patches = plt.hist(
        reaction_fluxes, bins=n_bins, density=False, facecolor="red", ec="black"
    )

    plt.xlabel("Flux (mmol/gDW/h)", fontsize=16)
    plt.ylabel("Frequency (#samples: " + str(reaction_fluxes.size) + ")", fontsize=14)
    plt.grid(True)
    plt.title("Reaction: " + reaction, fontweight="bold", fontsize=18)
    plt.axis([np.amin(reaction_fluxes), np.amax(reaction_fluxes), 0, np.amax(n) * 1.2])

    plt.show()
