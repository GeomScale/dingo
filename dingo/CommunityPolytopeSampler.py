# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file


import numpy as np
import math
from dingo.fva import slow_fva
from dingo.utils import (
    map_samples_to_steady_states,
    get_matrices_of_low_dim_polytope,
    get_matrices_of_full_dim_polytope,
    buildConqMatrix
)

try:
    import gurobipy
    from dingo.gurobi_based_implementations import fast_fba, fast_fva, fast_inner_ball
except ImportError as e:
    pass

from volestipy import HPolytope



class CommunityPolytopeSampler:

    def __init__(self, metabol_net):

        self._metabolic_network = metabol_net
        self._modelList          = metabol_net.modelList
        self.list_of_model_elements = []
        self._comm_A       = []
        self._comm_b       = []
        self._comm_N       = []
        self._comm_N_shift = []
        self._comm_T       = []
        self._comm_T_shift = []

        self._parameters = {}
        self._parameters["nullspace_method"] = "sparseQR"
        self._parameters["opt_percentage"] = self.metabolic_network.parameters[
            "opt_percentage"
        ]
        self._parameters["distribution"] = "uniform"
        self._parameters["first_run_of_mmcs"] = True

        try:
            import gurobipy

            self._parameters["fast_computations"] = True
            self._parameters["tol"] = 1e-06
        except ImportError as e:
            self._parameters["fast_computations"] = False
            self._parameters["tol"] = 1e-03


    def getIndividualMatrices(self):
        """
        A Python function to derive the matrices A, Aeq and the vectors b, beq for each model. 
        Here is what each of these variables stand for:

        A   -- Linear inequality constraints, specified as a real matrix. A is an M-by-N matrix, where M is the number of inequalities, and N is the number of variables 
        b   -- Linear inequality constraints, specified as a real vector. b is an M-element vector related to the A matrix. I
        Aeq -- Linear equality constraints, specified as a real matrix. Aeq is an Me-by-N matrix, where Me is the number of equalities, and N is the number of variables
        beq -- Linear equality constraints, specified as a real vector. beq is an Me-element vector related to the Aeq matrix.
        min_fluxes
        max_fluxes
        """
        list_of_model_elements = []
        for model in self._modelList:

            polytope_cl = PolytopeSampler(model)

            if (
                polytope_cl.A == []
                or polytope_cl.b == []
                or polytope_cl.N == []
                or polytope_cl.N_shift == []
                or polytope_cl.T == []
                or polytope_cl.T_shift == []
            ):

                (
                    min_fluxes,
                    max_fluxes,
                    max_biomass_flux_vector,
                    max_biomass_objective,
                ) = polytope_cl.metabolic_network.fva()

                A, b, Aeq, beq = get_matrices_of_low_dim_polytope(
                    polytope_cl.metabolic_network.S,
                    polytope_cl.metabolic_network.lb,
                    polytope_cl.metabolic_network.ub,
                    min_fluxes,
                    max_fluxes,
                )

            # Make a tupple with all needed for each model
            model_elements = (A, b, Aeq, beq, min_fluxes, max_fluxes)
            list_of_model_elements.append(model_elements)

        self.list_of_model_elements = list_of_model_elements

    def matrices_for_community_level(self):
        """A Python function to derive the matrices A, Aeq and the vectors b, beq of the low dimensional polytope at the community level,
        such that A*x <= b and Aeq*x = beq.

        Keyword arguments:
        self.list_of_model_elements -- output of the getIndividualMatrices function, including
        A, b, Aeq, beq, min_fluxes, max_fluxes for each model 
        """

        list_of_A   = []
        list_of_b   = []
        list_of_Aeq = []
        list_of_beq = []

        for model in self.list_of_model_elements:

            list_of_A.append(model[0])
            list_of_b.append(model[1])
            list_of_Aeq.append(model[2])
            list_of_beq.append (model[3])

        tmp_A   = buildConqMatrix(list_of_A)
        tmp_b   = np.concatenate(list_of_b, axis=0)

        tmp_Aeq = buildConqMatrix(list_of_Aeq)
        tmp_beq = np.concatenate(list_of_beq, axis=0)


        # By making use of the matrices just built, get full polytope  
        (
            self._comm_A, 
            self._comm_b, 
            self._comm_N, 
            self._comm_N_shift, 
        ) = get_matrices_of_full_dim_polytope(tmp_A, tmp_b, tmp_Aeq, tmp_beq)

        n = self._comm_A.shape[1]
        self._T = np.eye(n)
        self._T_shift = np.zeros(n)

        return self._comm_A, self._comm_b, self._comm_N, self._comm_N_shift


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

        self.getIndividualMatrices()
        self.matrices_for_community_level()

        P = HPolytope(self._comm_A, self._comm_b)

        if self._parameters["fast_computations"]:
            self._comm_A, self._comm_b, Tr, Tr_shift, samples = P.fast_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )
        else:
            self._comm_A, self._comm_b, Tr, Tr_shift, samples = P.slow_mmcs(
                ess, psrf, parallel_mmcs, num_threads
            )

        if self._parameters["first_run_of_mmcs"]:
            steady_states = map_samples_to_steady_states(
                samples, self._comm_N, self._comm_N_shift
            )
            self._parameters["first_run_of_mmcs"] = False
        else:
            steady_states = map_samples_to_steady_states(
                samples, self._comm_N, self._comm_N_shift, self._T, self._T_shift
            )

        self._T = np.dot(self._T, Tr)
        self._T_shift = np.add(self._T_shift, Tr_shift)

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
