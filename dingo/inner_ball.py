# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp


# This is a function to get the maximum ball inscribed in the polytope using scipy.optimize LP solver `linprog`
def slow_inner_ball(A, b):
    """A Python function to compute the maximum inscribed ball in the given polytope using scipy.optimize LP solver `linprog`
    Returns the optimal solution for the following linear program:
    max r, subject to,
    a_ix + r||a_i|| <= b, i=1,...,n

    Keyword arguments:
    A -- an mxn matrix that contains the normal vectors of the facets of the polytope row-wise
    b -- a m-dimensional vector
    """

    extra_column = []

    m = A.shape[0]
    n = A.shape[1]

    for i in range(A.shape[0]):
        entry = np.linalg.norm(
            A[
                i,
            ]
        )
        extra_column.append(entry)

    column = np.asarray(extra_column)
    A_expand = np.c_[A, column]

    a = -np.ones((1, n + 1))
    azero = np.zeros((1, n))
    a[:, :-1] = azero
    objective_function = a[0]

    res = linprog(objective_function, A_ub=A_expand, b_ub=b, bounds=(None, None))

    sol = res.x

    # Get the center point and the radius of max ball from the solution of LP; its last element is the radius
    point = []
    for i in range(len(sol)):
        if i == len(sol) - 1:
            r = sol[i]
        else:
            value = sol[i]
            point.append(value)

    # And check whether the computed radius is negative
    if r < 0:
        raise Exception(
            "The full dimensional polytope is infeasible or the solver failed to compute the maximum inscribed ball."
        )
    else:
        return point, r
