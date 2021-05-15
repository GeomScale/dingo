# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Vissarion Fisikopoulos

# Licensed under GNU LGPL.3, see LICENCE file

import unittest
import os
import scipy
import numpy as np
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.inner_ball import slow_inner_ball
from dingo.scaling import gmscale


class TestMaxBall(unittest.TestCase):
    
    def test_simple(self):
        m = 2
        n = 5
        A = np.zeros((2 * n, n), dtype="float")
        A[0:n] = np.eye(n)
        A[n:] -= np.eye(n, n, dtype="float")
        b = np.ones(2 * n, dtype="float")

        max_ball = slow_inner_ball(A, b)

        self.assertTrue(abs(max_ball[1] - 1) < 1e-04)

if __name__ == "__main__":
    unittest.main()
