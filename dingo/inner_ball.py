# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
from scipy.optimize import linprog
import scipy.sparse as sp


# This is a function to get the maximum ball inscribed in the polytope using scipy.optimize LP solver `linprog`
def slow_inner_ball(A, b):

    extra_column = []

    m = A.shape[0]
    n = A.shape[1]

    for i in range(A.shape[0]):
        entry = np.linalg.norm(A[i,])
        extra_column.append(entry)

    column = np.asarray(extra_column)
    A_expand = np.c_[A, column]
    
    a = -np.ones((1,n+1))
    azero = np.zeros((1,n))
    a[:,:-1] = azero
    objective_function = a[0]

    res = linprog(objective_function, A_ub = A_expand, b_ub = b, bounds = (None, None))

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
    if r < 0 :
        print ("The radius calculated has negative value. The polytope is infeasible or something went wrong with the solver")
    else:
        return point, r

