# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
from scipy import linalg
import sparseqr
import scipy.sparse.linalg


# Build a Python function to compute the nullspace of the stoichiometric matrix and a shifting to the origin
def nullspace_dense(Aeq, beq):
    """A Python function to compute the matrix of the right nullspace of the augmented stoichiometric
    matrix and a shifting to the origin
    (a) Solves the equation Aeq x = beq,
    (b) Computes the nullspace of Aeq

    Keyword arguments:
    Aeq -- the mxn augmented row-wise stoichiometric matrix
    beq -- a m-dimensional vector
    """

    N_shift = np.linalg.lstsq(Aeq, beq, rcond=None)[0]
    N = linalg.null_space(Aeq)

    return N, N_shift


def nullspace_sparse(Aeq, beq):
    """A Python function to compute the matrix of the right nullspace of the augmented stoichiometric
    matrix, exploiting that the matrix is in sparse format, and a shifting to the origin.
    The function uses the python wrapper PySPQR for the SuiteSparseQR library to compute the QR decomposition of matrix Aeq
    (a) Solves the equation Aeq x = beq,
    (b) Computes the nullspace of Aeq

    Keyword arguments:
    Aeq -- the mxn augmented row-wise stoichiometric matrix
    beq -- a m-dimensional vector
    """

    N_shift = np.linalg.lstsq(Aeq, beq, rcond=None)[0]
    Aeq = Aeq.T
    Aeq = scipy.sparse.csc_matrix(Aeq)

    # compute the QR decomposition of the Aeq_transposed
    Q, R, E, rank = sparseqr.qr(Aeq)

    # convert the matrices to dense format
    Q = Q.todense()
    R = R.todense()
    Aeq = Aeq.todense()

    if rank == 0:

        # Take the first n columns of Q where n is the number of columns of Aeq
        N = Q[:, : Aeq.shape[1]]

    else:

        # Take the last n-r columns of Q to derive the right nullspace of Aeq
        N = Q[:, rank:]

    N = np.asarray(N, dtype="float")
    N = np.ascontiguousarray(N, dtype="float")

    return N, N_shift
