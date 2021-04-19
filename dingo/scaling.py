# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import scipy.sparse as sp
from scipy.sparse import diags
import math


def gmscale(A, scltol):
    """This function is a python translation of the matlab cobra script you may find here:
    https://github.com/opencobra/cobratoolbox/blob/master/src/analysis/subspaces/gmscale.m
    Computes a scaling qA for the matrix A such that the computations in the polytope P = {x | qAx <= b}
    are numerically more stable

    Keyword arguments:
    A -- a mxn matrix
    scltol -- should be in the range (0.0, 1.0) to declare the desired accuracy. The maximum accuracy corresponds to 1.
    """

    m = A.shape[0]
    n = A.shape[1]
    A = np.abs(A)
    A_t = A.T  ## We will work with the transpose matrix as numpy works on row based
    maxpass = 10
    aratio = 1e50
    damp = 1e-4
    small = 1e-8
    rscale = np.ones((m, 1))
    cscale = np.ones((n, 1))

    # Main loop
    for npass in range(maxpass):

        rscale[rscale == 0] = 1
        sparse = sp.csr_matrix(1.0 / rscale)
        r_diagonal_elements = np.asarray(sparse.todense())
        Rinv = diags(r_diagonal_elements[0].T, shape=(m, m)).toarray()

        SA = np.dot(Rinv, A)
        SA_T = SA.T

        J, I = SA_T.nonzero()
        V = SA.T[SA.T > 0]
        invSA = sp.csr_matrix((1.0 / V, (J, I)), shape=(n, m)).T

        cmax = np.max(SA, axis=0)
        cmin = np.max(invSA, axis=0).data
        cmin = 1.0 / (cmin + 2.2204e-16)

        sratio = np.max(cmax / cmin)

        if npass > 0:
            c_product = np.multiply(
                np.max((np.array((cmin.data, damp * cmax))), axis=0), cmax
            )
            cscale = np.sqrt(c_product)

        check = aratio * scltol

        if npass >= 2 and sratio >= check:
            break

        if npass == maxpass:
            break

        aratio = sratio

        # Set new row scales for the next pass.
        cscale[cscale == 0] = 1
        sparse = sp.csr_matrix(1.0 / cscale)
        c_diagonal_elements = np.asarray(sparse.todense())

        Cinv = diags(c_diagonal_elements[0].T, shape=(n, n)).toarray()
        SA = np.dot(A, Cinv)
        SA_T = SA.T

        J, I = SA_T.nonzero()
        V = SA.T[SA.T > 0]
        invSA = sp.csr_matrix((1.0 / V, (J, I)), shape=(n, m)).T

        rmax = np.max(SA, axis=1)
        rmin = np.max(invSA, axis=1).data
        tmp = rmin + 2.2204e-16
        rmin = 1.0 / tmp

        r_product = np.multiply(
            np.max((np.array((rmin.data, damp * rmax))), axis=0), rmax
        )
        rscale = np.sqrt(r_product)

    # End of main loop

    # Reset column scales so the biggest element in each scaled column will be 1.
    # Again, allow for empty rows and columns.
    rscale[rscale == 0] = 1
    sparse = sp.csr_matrix(1.0 / rscale)
    r_diagonal_elements = np.asarray(sparse.todense())
    Rinv = diags(r_diagonal_elements[0].T, shape=(m, m)).toarray()

    SA = np.dot(Rinv, A)
    SA_T = SA.T
    J, I = SA_T.nonzero()
    V = SA.T[SA.T > 0]

    cscale = np.max(SA, axis=0)
    cscale[cscale == 0] = 1

    return cscale, rscale
