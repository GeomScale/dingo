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


def apply_scaling(A, b, cs, rs):
    """A Python function to apply the scaling computed by the function `gmscale` to a convex polytope

    Keyword arguments:
    A -- an mxn matrix that contains the normal vectors of the facets of the polytope row-wise
    b -- a m-dimensional vector
    cs -- a scaling vector for the matrix A
    rs -- a scaling vector for the vector b
    """

    m = rs.shape[0]
    n = cs.shape[0]
    r_diagonal_matrix = diags(1 / rs, shape=(m, m)).toarray()
    c_diagonal_matrix = diags(1 / cs, shape=(n, n)).toarray()

    new_A = np.dot(r_diagonal_matrix, np.dot(A, c_diagonal_matrix))
    new_b = np.dot(r_diagonal_matrix, b)

    return new_A, new_b, c_diagonal_matrix


def remove_almost_redundant_facets(A, b):
    """A Python function to remove the facets of a polytope with norm smaller than 1e-06

    Keyword arguments:
    A -- an mxn matrix that contains the normal vectors of the facets of the polytope row-wise
    b -- a m-dimensional vector
    """

    new_A = []
    new_b = []

    for i in range(A.shape[0]):
        entry = np.linalg.norm(
            A[
                i,
            ]
        )
        if entry < 1e-06:
            continue
        else:
            new_A.append(A[i, :])
            new_b.append(b[i])

    new_A = np.array(new_A)
    new_b = np.array(new_b)

    return new_A, new_b


# Map the points samples on the (rounded) full dimensional polytope, back to the initial one to obtain the steady states of the metabolic network
def map_samples_to_steady_states(samples, N, N_shift, T=None, T_shift=None):
    """A Python function to map back to the initial space the sampled points from a full dimensional polytope derived by two
    linear transformation of a low dimensional polytope, to obtain the steady states of the metabolic network

    Keyword arguments:
    samples -- an nxN matrix that contains sample points column-wise
    N, N_shift -- the matrix and the vector of the linear transformation applied on the low dimensional polytope to derive the full dimensional polytope
    T, T_shift -- the matrix and the vector of the linear transformation applied on the full dimensional polytope
    """

    extra_2 = np.full((samples.shape[1], N.shape[0]), N_shift)
    if T is None or T_shift is None:
        steady_states = N.dot(samples) + extra_2.T
    else:
        extra_1 = np.full((samples.shape[1], samples.shape[0]), T_shift)
        steady_states = N.dot(T.dot(samples) + extra_1.T) + extra_2.T

    return steady_states
