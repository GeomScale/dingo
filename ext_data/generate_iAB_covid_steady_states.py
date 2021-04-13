# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import pickle
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

from dingo import HPolytope


current_directory = os.getcwd()
input_file_mat = current_directory + "/iAB_AMO1410_SARS_CoV_2.mat"

iAB_AMO1410_SARS_CoV_2 = read_mat_file(input_file_mat)

lb = iAB_AMO1410_SARS_CoV_2[0]
ub = iAB_AMO1410_SARS_CoV_2[1]
S = iAB_AMO1410_SARS_CoV_2[2]
metabolites = iAB_AMO1410_SARS_CoV_2[3]
reactions = iAB_AMO1410_SARS_CoV_2[4]

covid_biomass_index = iAB_AMO1410_SARS_CoV_2[5]
human_biomass_index = covid_biomass_index - 1
biomass_function = np.zeros(S.shape[1])

step = 0.1
human_coeff = 1
covid_coeff = 0

name = "iAB_AMO1410_SARS_CoV_2"


for i in range(12):

    if i == 11:
        biomass_function = np.zeros(S.shape[1])
    else:
        biomass_function[human_biomass_index] = human_coeff
        biomass_function[covid_biomass_index] = covid_coeff

    print("covid: " + str(human_coeff))
    print("human: " + str(covid_coeff))

    print("computing FVA...")
    fva_res = fast_fva(lb, ub, S, biomass_function)
    print("FVA completed")

    A = fva_res[0]
    b = fva_res[1]
    Aeq = fva_res[2]
    beq = fva_res[3]
    min_fluxes = fva_res[4]
    max_fluxes = fva_res[5]

    nullspace_res = nullspace_sparse(Aeq, beq)

    N = nullspace_res[0]
    N_shift = nullspace_res[1]

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

    A, b, T, T_shift, samples = p.fast_mmcs(1000, True)

    steady_states = map_samples_to_steady_states(samples, N, N_shift)

    polytope_info = (
        A,
        b,
        N,
        N_shift,
        T,
        T_shift,
        samples,
    )
    network_info = (
        min_fluxes,
        max_fluxes,
    )
    metabol_reaction = (
        metabolites,
        reactions,
    )
    steadystates = (steady_states,)

    polytope_info = polytope_info + (name,)

    with open(
        "dingo_rounded_polytope_"
        + name
        + "_human_"
        + str(human_coeff)
        + "_covid_"
        + str(covid_coeff),
        "wb",
    ) as dingo_polytope_file:
        pickle.dump(polytope_info, dingo_polytope_file)

    with open(
        "dingo_minmax_fluxes_"
        + name
        + "_human_"
        + str(human_coeff)
        + "_covid_"
        + str(covid_coeff),
        "wb",
    ) as dingo_network_file:
        pickle.dump(network_info, dingo_network_file)

    with open(
        "dingo_metabolites_reactions_"
        + name
        + "_human_"
        + str(human_coeff)
        + "_covid_"
        + str(covid_coeff),
        "wb",
    ) as dingo_metabolreactions_file:
        pickle.dump(metabol_reaction, dingo_metabolreactions_file)

    with open(
        "dingo_steady_states_"
        + name
        + "_human_"
        + str(human_coeff)
        + "_covid_"
        + str(covid_coeff),
        "wb",
    ) as dingo_steadystates_file:
        pickle.dump(steadystates, dingo_steadystates_file)

    human_coeff = human_coeff - step
    covid_coeff = covid_coeff + step
