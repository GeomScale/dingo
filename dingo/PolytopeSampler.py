# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file


import numpy as np
import math
from dingo.MetabolicNetwork import MetabolicNetwork
from dingo.fva import slow_fva
from dingo.utils import (
    map_samples_to_steady_states,
    get_matrices_of_low_dim_polytope,
    get_matrices_of_full_dim_polytope,
)

try:
    import gurobipy
    from dingo.gurobi_based_implementations import fast_fba, fast_fva, fast_inner_ball
except ImportError as e:
    pass

from volestipy import HPolytope


class PolytopeSampler:
    def __init__(self, metabol_net):

        if not isinstance(metabol_net, MetabolicNetwork):
            raise Exception("An unknown input object given for initialization.")

        self.metabolic_network = metabol_net
        self.A = []
        self.b = []
        self.N = []
        self.N_shift = []
        self.T = []
        self.T_shift = []
        self.parameters = {}
        self.parameters["nullspace_method"] = "sparseQR"
        self.parameters["opt_percentage"] = self.metabolic_network.parameters[
            "opt_percentage"
        ]
        self.parameters["distribution"] = "uniform"
        self.parameters["first_run_of_mmcs"] = True

        try:
            import gurobipy

            self.parameters["fast_computations"] = True
            self.parameters["tol"] = 1e-06
        except ImportError as e:
            self.parameters["fast_computations"] = False
            self.parameters["tol"] = 1e-03

    def get_polytope(self):
        """A member function to derive the corresponding full dimensional polytope
        and a isometric linear transformation.
        """

        if (
            self.A == []
            or self.b == []
            or self.N == []
            or self.N_shift == []
            or self.T == []
            or self.T_shift == []
        ):

            (
                min_fluxes,
                max_fluxes,
                max_biomass_flux_vector,
                max_biomass_objective,
            ) = self.metabolic_network.fva()

            A, b, Aeq, beq = get_matrices_of_low_dim_polytope(
                self.metabolic_network.S,
                self.metabolic_network.lb,
                self.metabolic_network.ub,
                min_fluxes,
                max_fluxes,
            )

            if (
                A.shape[0] != b.size
                or A.shape[1] != Aeq.shape[1]
                or Aeq.shape[0] != beq.size
            ):
                raise Exception("FVA failed.")

            A = np.vstack((A, -self.metabolic_network.biomass_function))

            b = np.append(
                b,
                -(self.parameters["opt_percentage"] / 100)
                * self.parameters["tol"]
                * math.floor(max_biomass_objective / self.parameters["tol"]),
            )

            self.A, self.b, self.N, self.N_shift = get_matrices_of_full_dim_polytope(
                A, b, Aeq, beq
            )

            n = self.A.shape[1]
            self.T = np.eye(n)
            self.T_shift = np.zeros(n)

        return self.A, self.b, self.N, self.N_shift

    def generate_steady_states(
        self, ess=1000, psrf=False, parallel_mmcs=False, num_threads=1
    ):
        """A member function to sample steady states.

        Keyword arguments:
        ess -- the target effective sample size
        psrf -- a boolean flag to request PSRF smaller than 1.1 for all marginal fluxes
        parallel_mmcs -- a boolean flag to request the parallel mmcs
        num_threads -- the number of threads to use for parallel mmcs
        """

        self.get_polytope()

        P = HPolytope(self.A, self.b)

        if self.parameters["fast_computations"]:
            self.A, self.b, Tr, Tr_shift, samples = P.fast_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )
        else:
            self.A, self.b, Tr, Tr_shift, samples = P.slow_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )

        if self.parameters["first_run_of_mmcs"]:
            steady_states = map_samples_to_steady_states(samples, self.N, self.N_shift)
            self.parameters["first_run_of_mmcs"] = False
        else:
            steady_states = map_samples_to_steady_states(
                samples, self.N, self.N_shift, self.T, self.T_shift
            )

        self.T = np.dot(self.T, Tr)
        self.T_shift = np.add(self.T_shift, Tr_shift)

        return steady_states

    @staticmethod
    def sample_from_polytope(
        A, b, ess=1000, psrf=False, parallel_mmcs=False, num_threads=1
    ):
        """A static function to sample from a full dimensional polytope.

        Keyword arguments:
        A -- an mxn matrix that contains the normal vectors of the facets of the polytope row-wise
        b -- a m-dimensional vector, s.t. A*x <= b
        ess -- the target effective sample size
        psrf -- a boolean flag to request PSRF smaller than 1.1 for all marginal fluxes
        parallel_mmcs -- a boolean flag to request the parallel mmcs
        num_threads -- the number of threads to use for parallel mmcs
        """

        P = HPolytope(A, b)

        if self.parameters["fast_computations"]:
            A, b, Tr, Tr_shift, samples = P.fast_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )
        else:
            A, b, Tr, Tr_shift, samples = P.slow_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )

        return samples

    @staticmethod
    def sample_from_fva_output(
        min_fluxes,
        max_fluxes,
        biomass_function,
        max_biomass_objective,
        S,
        opt_percentage=100,
        ess=1000,
        psrf=False,
        parallel_mmcs=False,
        num_threads=1,
    ):
        """A static function to sample steady states when the output of FVA is given.

        Keyword arguments:
        min_fluxes -- minimum values of the fluxes, i.e., a n-dimensional vector
        max_fluxes -- maximum values for the fluxes, i.e., a n-dimensional vector
        biomass_function -- the biomass objective function
        max_biomass_objective -- the maximum value of the biomass objective function
        S -- stoichiometric matrix
        opt_percentage -- consider solutions that give you at least a certain
                      percentage of the optimal solution (default is to consider
                      optimal solutions only)
        ess -- the target effective sample size
        psrf -- a boolean flag to request PSRF smaller than 1.1 for all marginal fluxes
        parallel_mmcs -- a boolean flag to request the parallel mmcs
        num_threads -- the number of threads to use for parallel mmcs
        """

        A, b, Aeq, beq = get_matrices_of_low_dim_polytope(
            S, min_fluxes, max_fluxes, opt_percentage, tol
        )

        A = np.vstack((A, -biomass_function))
        b = np.append(
            b,
            -(opt_percentage / 100)
            * self.parameters["tol"]
            * math.floor(max_biomass_objective / self.parameters["tol"]),
        )

        A, b, N, N_shift = get_matrices_of_full_dim_polytope(A, b, Aeq, beq)

        P = HPolytope(A, b)

        if self.parameters["fast_computations"]:
            A, b, Tr, Tr_shift, samples = P.fast_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )
        else:
            A, b, Tr, Tr_shift, samples = P.slow_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )

        steady_states = map_samples_to_steady_states(samples, N, N_shift)

        return steady_states

    def set_fast_mode(self):

        self.parameters["fast_computations"] = True
        self.parameters["tol"] = 1e-06

    def set_slow_mode(self):

        self.parameters["fast_computations"] = False
        self.parameters["tol"] = 1e-03

    def set_distribution(self, value):

        self.parameters["distribution"] = value

    def set_nullspace_method(self, value):

        self.parameters["nullspace_method"] = value

    def set_tol(self, value):

        self.parameters["tol"] = value

    def set_opt_percentage(self, value):

        self.parameters["opt_percentage"] = value
