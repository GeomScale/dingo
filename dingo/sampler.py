# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file


import numpy as np
from dingo.metabolic_network import metabolic_network
from dingo.fva import slow_fva
from dingo.utils import (
    apply_scaling,
    remove_almost_redundant_facets,
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

class network_sampler:

    def __init__(self, metabol_net):

        if isinstance(metabol_net, metabolic_network):

            self.metabolic_network = metabolic_network
            self.A = []
            self.b = []
            self.N = []
            self.N_shift = []
            self.T =[]
            self.T_shift = []
            self.parameters = {}
            self.parameters['nullspace_method'] = 'sparseQR'
            self.parameters['opt_percentage'] = 100
            self.parameters['distribution'] = 'uniform'
            self.parameters['fast_computations'] = False
            self.parameters['tol'] = 1e-03

        else:

            raise Exception("An unknown input object given for initialization.")
    
    
    def get_polytope(self):

        if (self.A == [] or self.b == [] or self.N == [] or self.N_shift = [] or self.T = [] or self.T_shift = []):
            
            m = S.shape[0]
            n = S.shape[1]

            self.T =np.zeros((2 * n, n), dtype="float")
            self.T[0:n] = np.eye(n)
            self.T[n:] -= np.eye(n, n, dtype="float")

            self.T_shift = np.zeros(n)

            fva_res = fast_fva(self.metabolic_network.lb, self.metabolic_network.ub, self.metabolic_network.S, self.metabolic_network.biomass_function)

            AA = fva_res[0] # for testing
            bb = fva_res[1] # for testing
            Aeqq = fva_res[2] # for testing
            beqq = fva_res[3] # for testing
            min_fluxes = fva_res[4]
            max_fluxes = fva_res[5]
            max_biomass_flux_vector = fva_res[6]
            max_biomass_objective = fva_res[7]
            del fva_res

            A, b, Aeq, beq = get_matrices_of_low_dim_polytope(self.metabolic_network.S, self.metabolic_network.min_fluxes, self.metabolic_network.max_fluxes, self.parameters['opt_percentage'], self.parameters['tol'])

            if A.shape[0] != b.size or A.shape[1] != Aeq.shape[1] or Aeq.shape[0] != beq.size:
                raise Exception("FVA failed.")

            A = np.vstack((A, -self.metabolic_network.biomass_function))

            b = np.append(
                b, -(self.parameters['opt_percentage'] / 100) * self.parameters['tol'] * math.floor(max_biomass_objective / self.parameters['tol'])
            )

            self.A, self.b, self.N, self.N_shift = get_matrices_of_full_dim_polytope(A, b, Aeq, beq)

        return self.A, self.b, self.N, self.N_shift


    def generate_steady_states(self, ess = 1000, psrf = False, parallel_mmcs = False, num_threads = 1):

        get_polytope(self)

        P = HPolytope(self.A, self.b)

        self.A, self.b, Tr, Tr_shift, samples = P.fast_mmcs(ess, psrf, parallel_mmcs, num_threads)
        steady_states = map_samples_to_steady_states(samples, self.N, self.N_shift, self.T, self.T_shift)

        self.T = np.dot(self.T, Tr)
        self.T_shift = np.add(self.T_shift, Tr_shift)

        return steady_states

    @staticmethod
    def sample_from_polytope(A, b, ess = 1000, psrf = False, parallel_mmcs = False, num_threads = 1):

        P = HPolytope(A, b)

        A, b, Tr, Tr_shift, samples = P.fast_mmcs(ess, psrf, parallel_mmcs, num_threads)

        return samples

    @staticmethod
    def sample_from_fva_output(min_fluxes, max_fluxes, biomass_function, max_biomass_objective, S, opt_percentage = 100, tol = 1e-06, ess = 1000, psrf = False, parallel_mmcs = False, num_threads = 1):

        A, b, Aeq, beq = get_matrices_of_low_dim_polytope(S, min_fluxes, max_fluxes, opt_percentage, tol)

        A = np.vstack((A, -biomass_function))
        b = np.append(
            b, -(opt_percentage / 100) * tol * math.floor(max_biomass_objective / tol)
        )
        
        A, b, N, N_shift = get_matrices_of_full_dim_polytope(A, b, Aeq, beq)

        P = HPolytope(A, b)

        A, b, Tr, Tr_shift, samples = P.fast_mmcs(ess, psrf, parallel_mmcs, num_threads)
        steady_states = map_samples_to_steady_states(samples, N, N_shift)

        return steady_states
    
    def set_fast_mode(self):

        self.parameters['fast_computations'] = True
        self.parameters['tol'] = 1e-06

    def set_distribution(self, value):

        self.parameters['distribution'] = value
    
    def set_nullspace_method(self, value):

        self.parameters['nullspace_method'] = value

    def set_tol(self, value):

        self.parameters['tol'] = value

    def set_opt_percentage(self, value):

        self.parameters['opt_percentage'] = value
