# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file


import copy
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
import pandas as pd

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

        self._metabolic_network = copy.deepcopy(metabol_net)
        self._A = np.empty( shape=(0, 0) )
        self._b = np.empty( shape=(0, 0) )
        self._N = np.empty( shape=(0, 0) )
        self._N_shift = np.empty( shape=(0, 0) )
        self._T = np.empty( shape=(0, 0) )
        self._T_shift = np.empty( shape=(0, 0) )
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
            self._A.size == 0
            or self._b.size == 0
            or self._N.size == 0
            or self._N_shift.size == 0
            or self._T.size == 0
            or self._T_shift.size == 0
        ):

            (
                max_biomass_flux_vector,
                max_biomass_objective,
            ) = self._metabolic_network._fba()

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
                ) = self._metabolic_network._fva()

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

        steady_states_df = pd.DataFrame(steady_states, index = self._metabolic_network.reactions)

        return steady_states_df

    def generate_steady_states_no_multiphase(
        self, method = 'billiard_walk', n=1000, burn_in=0, thinning=1, variance=1.0, bias_vector=None
    ):
        """A member function to sample steady states.

        Keyword arguments:
        method -- An MCMC method to sample, i.e. {'billiard_walk', 'cdhr', 'rdhr', 'ball_walk', 'dikin_walk', 'john_walk', 'vaidya_walk', 'gaussian_hmc_walk', 'exponential_hmc_walk', 'hmc_leapfrog_gaussian', 'hmc_leapfrog_exponential'}
        n -- the number of steady states to sample
        burn_in -- the number of points to burn before sampling
        thinning -- the walk length of the chain
        """

        self.get_polytope()

        P = HPolytope(self._A, self._b)

        if bias_vector is None:
            bias_vector = np.ones(self._A.shape[1], dtype=np.float64)
        else:
            bias_vector = bias_vector.astype('float64')

        samples = P.generate_samples(method, n, burn_in, thinning, self._parameters["fast_computations"], variance, bias_vector)
        samples_T = samples.T

        steady_states = map_samples_to_steady_states(
                samples_T, self._N, self._N_shift
            )
        steady_states_df = pd.DataFrame(steady_states, index = self._metabolic_network.reactions)

        return steady_states_df

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
    def sample_from_polytope_no_multiphase(
        A, b, method = 'billiard_walk', n=1000, burn_in=0, thinning=1, variance=1.0, bias_vector=None
    ):
        """A static function to sample from a full dimensional polytope with an MCMC method.

        Keyword arguments:
        A -- an mxn matrix that contains the normal vectors of the facets of the polytope row-wise
        b -- a m-dimensional vector, s.t. A*x <= b
        method -- An MCMC method to sample, i.e. {'billiard_walk', 'cdhr', 'rdhr', 'ball_walk', 'dikin_walk', 'john_walk', 'vaidya_walk', 'gaussian_hmc_walk', 'exponential_hmc_walk', 'hmc_leapfrog_gaussian', 'hmc_leapfrog_exponential'}
        n -- the number of steady states to sample
        burn_in -- the number of points to burn before sampling
        thinning -- the walk length of the chain
        """
        if bias_vector is None:
            bias_vector = np.ones(A.shape[1], dtype=np.float64)
        else:
            bias_vector = bias_vector.astype('float64')

        P = HPolytope(A, b)

        try:
            import gurobipy
            samples = P.generate_samples(method, n, burn_in, thinning, True, variance, bias_vector)
        except ImportError as e:
            samples = P.generate_samples(method, n, burn_in, thinning, False, variance, bias_vector)

        samples_T = samples.T
        return samples_T

    @staticmethod
    def round_polytope(
        A, b, method = "john_position"
    ):
        P = HPolytope(A, b)
        try:
            import gurobipy
            A, b, Tr, Tr_shift, round_value = P.rounding(method, True)
        except ImportError as e:
            A, b, Tr, Tr_shift, round_value = P.rounding(method, False)

        return A, b, Tr, Tr_shift


    @staticmethod
    def sample_from_fva_output(
        model,
        opt_percentage=100,
        ess=1000,
        psrf=False,
        parallel_mmcs=False,
        num_threads=1,
    ):
        """A static function to sample steady states when the output of FVA is given.

        Keyword arguments:
        model -- a dingo.MetabolicNetwork() object
        opt_percentage -- consider solutions that give you at least a certain
                      percentage of the optimal solution (default is to consider
                      optimal solutions only)
        ess -- the target effective sample size
        psrf -- a boolean flag to request PSRF smaller than 1.1 for all marginal fluxes
        parallel_mmcs -- a boolean flag to request the parallel mmcs
        num_threads -- the number of threads to use for parallel mmcs
        """

        min_fluxes, max_fluxes, opt_vector, opt_value = model._fva()

        A, b, Aeq, beq = get_matrices_of_low_dim_polytope(
            model.S, min_fluxes, max_fluxes, opt_percentage, model._parameters["tol"]
        )

        A = np.vstack((A, -model.biomass_function))
        b = np.append(
            b,
            -(opt_percentage / 100)
            * model._parameters["tol"]
            * math.floor(opt_value / model._parameters["tol"]),
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

    @staticmethod
    def samples_as_df(model, samples):
        """A static function to convert the samples numpy ndarray to a pandas DataFrame with model's reactions as indices

        Keyword arguments:
        model --
        samples --
        """
        samples_df = pd.DataFrame(samples, index = model.reactions)
        return samples_df

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

