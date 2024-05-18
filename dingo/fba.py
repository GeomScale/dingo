# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp


def slow_fba(lb, ub, S, c):
    """A Python function to perform fba using scipy.optimize LP solver `linprog`
    Returns an optimal solution and its value for the following linear program:
    max c*v, subject to,
    Sv = 0, lb <= v <= ub

    Keyword arguments:
    lb -- lower bounds for the fluxes, i.e., a n-dimensional vector
    ub -- upper bounds for the fluxes, i.e., a n-dimensional vector
    S -- the mxn stoichiometric matrix, s.t. Sv = 0
    c -- the linear objective function, i.e., a n-dimensional vector
    """

    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )
    if c.size != S.shape[1]:
        raise Exception(
            "The length of the lineart objective function must be equal to the number of reactions."
        )

    m = S.shape[0]
    n = S.shape[1]
    optimum_value = 0
    optimum_sol = np.zeros(n)

    A = np.zeros((2 * n, n), dtype="float")
    A[0:n] = np.eye(n)
    A[n:] -= np.eye(n, n, dtype="float")

    b = np.concatenate((ub, -lb), axis=0)
    b = np.asarray(b, dtype="float")
    b = np.ascontiguousarray(b, dtype="float")

    beq = np.zeros(m)

    try:

        objective_function = np.asarray(c)
        objective_function = np.asarray([-x for x in c])

        res = linprog(
            objective_function, A_ub=A, b_ub=b, A_eq=S, b_eq=beq, bounds=(None, None)
        )

        # If optimized
        if res.success:
            optimum_value = -res.fun
            optimum_sol = np.asarray(res.x)

        return optimum_sol, optimum_value

    except AttributeError:
        print("scipy.optimize.linprog failed.")
