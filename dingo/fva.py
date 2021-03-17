# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
from scipy.optimize import linprog


def slow_fva(lb, ub, S):
    """A Python function to perform fva using scipy.optimize LP solver `linprog`
    Returns the value of the optimal solution for all the following linear programs:
    min/max v_i, for all coordinates i=1,...,n, subject to,
    Sv = 0, lb <= v <= ub

    Keyword arguments:
    lb -- lower bounds for the fluxes, i.e., a n-dimensional vector
    ub -- upper bounds for the fluxes, i.e., a n-dimensional vector
    S -- the mxn stoichiometric matrix, s.t. Sv = 0
    """

    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )

    m = S.shape[0]
    n = S.shape[1]

    A = np.zeros((2 * n, n), dtype="float")
    A[0:n] = np.eye(n)
    A[n:] -= np.eye(n, n, dtype="float")

    b = np.concatenate((ub, -lb), axis=0)
    b = np.asarray(b, dtype="float")
    b = np.ascontiguousarray(b, dtype="float")
    beq = np.zeros(m)

    Aeq_new = S
    beq_new = beq

    min_fluxes = []
    max_fluxes = []

    try:

        for i in range(int(A.shape[0] / 2)):

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

            # Calculate the width
            width = abs(max_objective - min_objective)

            # Check whether we keep or not the equality
            if width < 1e-07:
                Aeq_new = np.vstack(
                    (
                        Aeq_new,
                        A[
                            i,
                        ],
                    )
                )
                beq_new = np.append(beq_new, min(max_objective, min_objective))

        # The np.vstack() creates issues on changing contiguous c orded of np arrays; here we fix this
        Aeq_new = np.ascontiguousarray(Aeq_new, dtype=np.dtype)

        # Furthremore, we need to have float64 in all numpy arrays
        Aeq_new = Aeq_new.astype("float64")

        # Make lists of fluxes numpy arrays
        min_fluxes = np.asarray(min_fluxes)
        max_fluxes = np.asarray(max_fluxes)

        return A, b, Aeq_new, beq_new, min_fluxes, max_fluxes

    except AttributeError:
        print("scipy.optimize.linprog failed.")
