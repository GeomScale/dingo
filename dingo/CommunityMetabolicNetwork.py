# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import sys, os
from dingo.loading_models import read_json_file, read_mat_file, getModelList
from dingo.fva import slow_fva
from dingo.fba import slow_fba

try:
    import gurobipy
    from dingo.gurobi_based_implementations import fast_fba, fast_fva, fast_inner_ball
except ImportError as e:
    pass

class CommunityMetabolicNetwork():

    # This implementation works only for communities of 2 models 
    # Once our method is validated, we will move on to implement it for more 

    def __init__(self, comm_tuple_args):

        self._comm_parameters = {}
        self._comm_parameters["opt_percentage"]   = 100
        self._comm_parameters["distribution"]     = "uniform"
        self._comm_parameters["nullspace_method"] = "sparseQR"

        try:
            import gurobipy

            self._comm_parameters["fast_computations"] = True
        except ImportError as e:
            self._comm_parameters["fast_computations"] = False

        if len(comm_tuple_args) != 8:
            raise Exception(
                "An unknown input format given to initialize a metabolic network object."
            )

        self._comm_lb               = comm_tuple_args[0]
        self._comm_ub               = comm_tuple_args[1]
        self._comm_S                = comm_tuple_args[2]
        self._comm_metabolites      = comm_tuple_args[3]
        self._comm_reactions        = comm_tuple_args[4]
        self._comm_biomass_index    = comm_tuple_args[5]
        self._comm_biomass_function = comm_tuple_args[6]
        self._modelList        = comm_tuple_args[7]

        try:
            if (
                self._comm_lb.size                  != self._comm_ub.size
                or self._comm_lb.size               != self._comm_S.shape[1]
                or len(self._comm_metabolites)      != self._comm_S.shape[0]
                or len(self._comm_reactions)        != self._comm_S.shape[1]
                or self._comm_biomass_function.size != self._comm_S.shape[1]
                or (self._comm_biomass_index < 0)
                or (self._comm_biomass_index > self._comm_biomass_function.size)
            ):
                raise Exception(
                    "Wrong tuple format given to initialize a metabolic network object."
                )
        except LookupError as error:
            raise error.with_traceback(sys.exc_info()[2])


    @classmethod
    def buildModelList(cls, directory, format_type):
        comm_tuple_args = getModelList(directory, format_type)
        return cls(comm_tuple_args)

    def fva(self):
        """A member function to apply the FVA method on the community metabolic network."""

        if self._comm_parameters["fast_computations"]:
            return fast_fva(
                self._comm_lb,
                self._comm_ub,
                self._comm_S,
                self._comm_biomass_function,
                self._comm_parameters["opt_percentage"],
            )
        else:
            return slow_fva(
                self._comm_lb,
                self._comm_ub,
                self._comm_S,
                self._comm_biomass_function,
                self._comm_parameters["opt_percentage"],
            )

    def fba(self):
        """A member function to apply the FBA method on the community metabolic network."""

        if self._comm_parameters["fast_computations"]:
            return fast_fba(self._comm_lb, self._comm_ub, self._S, self._comm_biomass_function)
        else:
            return slow_fba(self._comm_lb, self._comm_ub, self._S, self._comm_biomass_function)

    @property
    def modelList(self):
        return self._modelList

    @property
    def lb(self):
        return self._comm_lb

    @property
    def ub(self):
        return self._comm_ub

    @property
    def S(self):
        return self._comm_S

    @property
    def metabolites(self):
        return self._comm_metabolites

    @property
    def reactions(self):
        return self._comm_reactions

    @property
    def biomass_index(self):
        return self._comm_biomass_index

    @property
    def biomass_function(self):
        return self._comm_biomass_function

    @property
    def parameters(self):
        return self._comm_parameters

    @property
    def get_as_tuple(self):
        return (
            self._comm_lb,
            self._comm_ub,
            self._comm_S,
            self._comm_metabolites,
            self._comm_reactions,
            self._comm_biomass_index,
            self._comm_biomass_function,
        )

    def num_of_reactions(self):
        return len(self._comm_reactions)

    def num_of_metabolites(self):
        return len(self._comm_metabolites)

    @lb.setter
    def lb(self, value):
        self._comm_lb = value

    @ub.setter
    def ub(self, value):
        self._comm_ub = value

    @S.setter
    def S(self, value):
        self._comm_S = value

    @metabolites.setter
    def metabolites(self, value):
        self._comm_metabolites = value

    @reactions.setter
    def reactions(self, value):
        self._comm_reactions = value

    @biomass_index.setter
    def biomass_index(self, value):
        self._comm_biomass_index = value

    @biomass_function.setter
    def biomass_function(self, value):
        self._comm_biomass_function = value

    def set_fast_mode(self):

        try:
            import gurobipy

            self._comm_parameters["fast_computations"] = True
        except ImportError as e:
            print("You have to install gurobi to use the fast computations.")
            self._comm_parameters["fast_computations"] = False

    def set_slow_mode(self):

        self._comm_parameters["fast_computations"] = False

    def set_nullspace_method(self, value):

        self._comm_parameters["nullspace_method"] = value

    def set_opt_percentage(self, value):

        self._comm_parameters["opt_percentage"] = value

    def shut_down_reaction(self, index_val):

        if (
            (not isinstance(index_val, int))
            or index_val < 0
            or index_val >= self._comm_S.shape[1]
        ):
            raise Exception("The input does not correspond to a proper reaction index.")

        self._comm_lb[index_val] = 0
        self._comm_ub[index_val] = 0 
