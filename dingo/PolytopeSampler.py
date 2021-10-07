# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file


import numpy as np
import warnings
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
    from dingo.gurobi_based_implementations import (
        fast_fba,
        fast_fva,
        fast_inner_ball,
        fast_remove_redundant_facets,
    )
except ImportError as e:
    pass

from volestipy import HPolytope


class PolytopeSampler:
    def __init__(self, metabol_net):

        if not isinstance(metabol_net, MetabolicNetwork):
            raise Exception("An unknown input object given for initialization.")

        self._metabolic_network = metabol_net
        self._A = []
        self._b = []
        self._N = []
        self._N_shift = []
        self._T = []
        self._T_shift = []
        self._parameters = {}
        self._parameters["nullspace_method"] = "sparseQR"
        self._parameters["opt_percentage"] = self.metabolic_network.parameters[
            "opt_percentage"
        ]
        self._parameters["distribution"] = "uniform"
        self._parameters["first_run_of_mmcs"] = True
        self._parameters["remove_redundant_facets"] = True

        try:
            import gurobipy

            self._parameters["fast_computations"] = True
            self._parameters["tol"] = 1e-06
        except ImportError as e:
            self._parameters["fast_computations"] = False
            self._parameters["tol"] = 1e-03

    def get_polytope(self):
        """A member function to derive the corresponding full dimensional polytope
        and a isometric linear transformation that maps the latter to the initial space.
        """

        if (
            self._A == []
            or self._b == []
            or self._N == []
            or self._N_shift == []
            or self._T == []
            or self._T_shift == []
        ):

            (
                max_biomass_flux_vector,
                max_biomass_objective,
            ) = self._metabolic_network.fba()

            if (
                self._parameters["fast_computations"]
                and self._parameters["remove_redundant_facets"]
            ):

                A, b, Aeq, beq = fast_remove_redundant_facets(
                    self._metabolic_network.lb,
                    self._metabolic_network.ub,
                    self._metabolic_network.S,
                    self._metabolic_network.biomass_function,
                    self._parameters["opt_percentage"],
                )
            else:
                if (not self._parameters["fast_computations"]) and self._parameters[
                    "remove_redundant_facets"
                ]:
                    warnings.warn(
                        "We continue without redundancy removal (slow mode is ON)"
                    )

                (
                    min_fluxes,
                    max_fluxes,
                    max_biomass_flux_vector,
                    max_biomass_objective,
                ) = self._metabolic_network.fva()

                A, b, Aeq, beq = get_matrices_of_low_dim_polytope(
                    self._metabolic_network.S,
                    self._metabolic_network.lb,
                    self._metabolic_network.ub,
                    min_fluxes,
                    max_fluxes,
                )

            if (
                A.shape[0] != b.size
                or A.shape[1] != Aeq.shape[1]
                or Aeq.shape[0] != beq.size
            ):
                raise Exception("Preprocess for full dimensional polytope failed.")

            A = np.vstack((A, -self._metabolic_network.biomass_function))

            b = np.append(
                b,
                -np.floor(max_biomass_objective / self._parameters["tol"])
                * self._parameters["tol"]
                * self._parameters["opt_percentage"]
                / 100,
            )

            (
                self._A,
                self._b,
                self._N,
                self._N_shift,
            ) = get_matrices_of_full_dim_polytope(A, b, Aeq, beq)

            n = self._A.shape[1]
            self._T = np.eye(n)
            self._T_shift = np.zeros(n)

        return self._A, self._b, self._N, self._N_shift

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

        P = HPolytope(self._A, self._b)

        if self._parameters["fast_computations"]:
            self._A, self._b, Tr, Tr_shift, samples = P.fast_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )
        else:
            self._A, self._b, Tr, Tr_shift, samples = P.slow_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )

        if self._parameters["first_run_of_mmcs"]:
            steady_states = map_samples_to_steady_states(
                samples, self._N, self._N_shift
            )
            self._parameters["first_run_of_mmcs"] = False
        else:
            steady_states = map_samples_to_steady_states(
                samples, self._N, self._N_shift, self._T, self._T_shift
            )

        self._T = np.dot(self._T, Tr)
        self._T_shift = np.add(self._T_shift, Tr_shift)

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

        try:
            import gurobipy

            A, b, Tr, Tr_shift, samples = P.fast_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )
        except ImportError as e:
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
            * self._parameters["tol"]
            * math.floor(max_biomass_objective / self._parameters["tol"]),
        )

        A, b, N, N_shift = get_matrices_of_full_dim_polytope(A, b, Aeq, beq)

        P = HPolytope(A, b)

        try:
            import gurobipy

            A, b, Tr, Tr_shift, samples = P.fast_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )
        except ImportError as e:
            A, b, Tr, Tr_shift, samples = P.slow_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )

        steady_states = map_samples_to_steady_states(samples, N, N_shift)

        return steady_states

    @property
    def A(self):
        return self._A

    @property
    def b(self):
        return self._b

    @property
    def T(self):
        return self._T

    @property
    def T_shift(self):
        return self._T_shift

    @property
    def N(self):
        return self._N

    @property
    def N_shift(self):
        return self._N_shift

    @property
    def metabolic_network(self):
        return self._metabolic_network

    def facet_redundancy_removal(self, value):
        self._parameters["remove_redundant_facets"] = value

        if (not self._parameters["fast_computations"]) and value:
            warnings.warn(
                "Since you are in slow mode the redundancy removal step is skipped (dingo does not currently support this functionality in slow mode)"
            )

    def set_fast_mode(self):

        self._parameters["fast_computations"] = True
        self._parameters["tol"] = 1e-06

    def set_slow_mode(self):

        self._parameters["fast_computations"] = False
        self._parameters["tol"] = 1e-03

    def set_distribution(self, value):

        self._parameters["distribution"] = value

    def set_nullspace_method(self, value):

        self._parameters["nullspace_method"] = value

    def set_tol(self, value):

        self._parameters["tol"] = value

    def set_opt_percentage(self, value):

        self._parameters["opt_percentage"] = value
