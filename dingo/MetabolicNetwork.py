# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Vissarion Fisikopoulos

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import sys
from dingo.loading_models import read_json_file, read_mat_file, read_sbml_file
from dingo.fva import slow_fva
from dingo.fba import slow_fba

try:
    import gurobipy
    from dingo.gurobi_based_implementations import fast_fba, fast_fva, fast_inner_ball
except ImportError as e:
    pass


class MetabolicNetwork:
    def __init__(self, tuple_args):

        self._parameters = {}
        self._parameters["opt_percentage"] = 100
        self._parameters["distribution"] = "uniform"
        self._parameters["nullspace_method"] = "sparseQR"

        try:
            import gurobipy

            self._parameters["fast_computations"] = True
        except ImportError as e:
            self._parameters["fast_computations"] = False

        if len(tuple_args) != 7:
            raise Exception(
                "An unknown input format given to initialize a metabolic network object."
            )

        self._lb = tuple_args[0]
        self._ub = tuple_args[1]
        self._S = tuple_args[2]
        self._metabolites = tuple_args[3]
        self._reactions = tuple_args[4]
        self._biomass_index = tuple_args[5]
        self._biomass_function = tuple_args[6]

        try:
            if self._biomass_index is not None and (
                self._lb.size != self._ub.size
                or self._lb.size != self._S.shape[1]
                or len(self._metabolites) != self._S.shape[0]
                or len(self._reactions) != self._S.shape[1]
                or self._biomass_function.size != self._S.shape[1]
                or (self._biomass_index < 0)
                or (self._biomass_index > self._biomass_function.size)
            ):
                raise Exception(
                    "Wrong tuple format given to initialize a metabolic network object."
                )
        except LookupError as error:
            raise error.with_traceback(sys.exc_info()[2])

    @classmethod
    def from_json(cls, arg):
        if (not isinstance(arg, str)) or (arg[-4:] != "json"):
            raise Exception(
                "An unknown input format given to initialize a metabolic network object."
            )

        return cls(read_json_file(arg))

    @classmethod
    def from_mat(cls, arg):
        if (not isinstance(arg, str)) or (arg[-3:] != "mat"):
            raise Exception(
                "An unknown input format given to initialize a metabolic network object."
            )

        return cls(read_mat_file(arg))

    @classmethod
    def from_sbml(cls, arg):
        if (not isinstance(arg, str)) or (arg[-3:] != "xml"):
            raise Exception(
                "An unknown input format given to initialize a metabolic network object."
            )

        return cls(read_sbml_file(arg))

    def fva(self):
        """A member function to apply the FVA method on the metabolic network."""

        if self._parameters["fast_computations"]:
            return fast_fva(
                self._lb,
                self._ub,
                self._S,
                self._biomass_function,
                self._parameters["opt_percentage"],
            )
        else:
            return slow_fva(
                self._lb,
                self._ub,
                self._S,
                self._biomass_function,
                self._parameters["opt_percentage"],
            )

    def fba(self):
        """A member function to apply the FBA method on the metabolic network."""

        if self._parameters["fast_computations"]:
            return fast_fba(self._lb, self._ub, self._S, self._biomass_function)
        else:
            return slow_fba(self._lb, self._ub, self._S, self._biomass_function)

    @property
    def lb(self):
        return self._lb

    @property
    def ub(self):
        return self._ub

    @property
    def S(self):
        return self._S

    @property
    def metabolites(self):
        return self._metabolites

    @property
    def reactions(self):
        return self._reactions

    @property
    def biomass_index(self):
        return self._biomass_index

    @property
    def biomass_function(self):
        return self._biomass_function

    @property
    def parameters(self):
        return self._parameters

    @property
    def get_as_tuple(self):
        return (
            self._lb,
            self._ub,
            self._S,
            self._metabolites,
            self._reactions,
            self._biomass_index,
            self._biomass_function,
        )

    def num_of_reactions(self):
        return len(self._reactions)

    def num_of_metabolites(self):
        return len(self._metabolites)

    @lb.setter
    def lb(self, value):
        self._lb = value

    @ub.setter
    def ub(self, value):
        self._ub = value

    @S.setter
    def S(self, value):
        self._S = value

    @metabolites.setter
    def metabolites(self, value):
        self._metabolites = value

    @reactions.setter
    def reactions(self, value):
        self._reactions = value

    @biomass_index.setter
    def biomass_index(self, value):
        self._biomass_index = value

    @biomass_function.setter
    def biomass_function(self, value):
        self._biomass_function = value

    def set_fast_mode(self):

        try:
            import gurobipy

            self._parameters["fast_computations"] = True
        except ImportError as e:
            print("You have to install gurobi to use the fast computations.")
            self._parameters["fast_computations"] = False

    def set_slow_mode(self):

        self._parameters["fast_computations"] = False

    def set_nullspace_method(self, value):

        self._parameters["nullspace_method"] = value

    def set_opt_percentage(self, value):

        self._parameters["opt_percentage"] = value

    def shut_down_reaction(self, index_val):

        if (
            (not isinstance(index_val, int))
            or index_val < 0
            or index_val >= self._S.shape[1]
        ):
            raise Exception("The input does not correspond to a proper reaction index.")

        self._lb[index_val] = 0
        self._ub[index_val] = 0
