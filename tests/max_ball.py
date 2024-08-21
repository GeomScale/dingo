# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Vissarion Fisikopoulos
# Copyright (c) 2024 Ke Shi

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import sys
import scipy
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.pyoptinterface_based_impl import inner_ball, set_default_solver
from dingo.scaling import gmscale


class TestMaxBall(unittest.TestCase):
    
    def test_simple(self):
        m = 2
        n = 5
        A = np.zeros((2 * n, n), dtype="float")
        A[0:n] = np.eye(n)
        A[n:] -= np.eye(n, n, dtype="float")
        b = np.ones(2 * n, dtype="float")

        max_ball = inner_ball(A, b)

        self.assertTrue(abs(max_ball[1] - 1) < 1e-04)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        set_default_solver(sys.argv[1])
        sys.argv.pop(1)
    unittest.main()
