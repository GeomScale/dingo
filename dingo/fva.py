# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
from scipy.optimize import linprog
from .fba import slow_fba
import math


def slow_fva(lb, ub, S, c, opt_percentage=100):
    """A Python function to perform fva using scipy.optimize LP solver `linprog`
    Returns the value of the optimal solution for all the following linear programs:
    min/max v_i, for all coordinates i=1,...,n, subject to,
    Sv = 0, lb <= v <= ub

    Keyword arguments:
    lb -- lower bounds for the fluxes, i.e., a n-dimensional vector
    ub -- upper bounds for the fluxes, i.e., a n-dimensional vector
    S -- the mxn stoichiometric matrix, s.t. Sv = 0
    c -- the objective function to maximize
    opt_percentage -- consider solutions that give you at least a certain
                      percentage of the optimal solution (default is to consider
                      optimal solutions only)
    """

    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )

    # declare the tolerance that linprog works properly (we found it experimentally)
    tol = 1e-03

    m = S.shape[0]
    n = S.shape[1]

    A = np.zeros((2 * n, n), dtype="float")
    A[0:n] = np.eye(n)
    A[n:] -= np.eye(n, n, dtype="float")

    b = np.concatenate((ub, -lb), axis=0)
    b = np.asarray(b, dtype="float")
    b = np.ascontiguousarray(b, dtype="float")
    beq = np.zeros(m)

    # call fba to obtain an optimal solution
    max_biomass_flux_vector, max_biomass_objective = slow_fba(lb, ub, S, c)

    # add an additional constraint to impose solutions with at least `opt_percentage` of the optimal solution
    A = np.vstack((A, -c))

    b = np.append(
        b, -(opt_percentage / 100) * tol * math.floor(max_biomass_objective / tol)
    )

    min_fluxes = []
    max_fluxes = []

    try:

        for i in range(n):

            # Set the ith row of the A matrix as the objective function
            objective_function = A[
                i,
            ]
            res = linprog(
                objective_function,
                A_ub=A,
                b_ub=b,
                A_eq=S,
                b_eq=beq,
                bounds=(None, None),
            )

            # If optimized
            if res.success:
                min_objective = res.fun
                min_fluxes.append(min_objective)
            else:
                min_fluxes.append(lb[i])

            # Set the minus of the ith row of the A matrix as the objective function
            objective_function = np.asarray([-x for x in objective_function])
            res = linprog(
                objective_function,
                A_ub=A,
                b_ub=b,
                A_eq=S,
                b_eq=beq,
                bounds=(None, None),
            )

            # If optimized
            if res.success:
                max_objective = -res.fun
                max_fluxes.append(max_objective)
            else:
                max_fluxes.append(ub[i])

        # Make lists of fluxes numpy arrays
        min_fluxes = np.asarray(min_fluxes)
        max_fluxes = np.asarray(max_fluxes)

        return (
            min_fluxes,
            max_fluxes,
            max_biomass_flux_vector,
            max_biomass_objective,
        )

    except AttributeError:
        print("scipy.optimize.linprog failed.")
