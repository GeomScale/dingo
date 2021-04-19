# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import sys
from dingo.loading_models import read_json_file, read_mat_file
from dingo.fva import slow_fva
from dingo.fba import slow_fba

try:
    import gurobipy
    from dingo.gurobi_based_implementations import fast_fba, fast_fva, fast_inner_ball
except ImportError as e:
    pass


class metabolic_network:
    def __init__(self, arg):

        self.parameters = {}
        self.parameters["opt_percentage"] = 100
        self.parameters["distribution"] = "uniform"
        self.parameters["nullspace_method"] = "sparseQR"

        try:
            import gurobipy

            self.parameters["fast_computations"] = True
        except ImportError as e:
            self.parameters["fast_computations"] = False

        if isinstance(arg, str):

            if arg[-3:] != "mat" and arg[-4:] != "json":
                raise Exception("An unknown format file given.")
            try:
                if arg[-4:] == "json":
                    (
                        self.lb,
                        self.ub,
                        self.S,
                        self.metabolites,
                        self.reactions,
                        self.biomass_index,
                        self.biomass_function,
                    ) = read_json_file(arg)
                else:
                    (
                        self.lb,
                        self.ub,
                        self.S,
                        self.metabolites,
                        self.reactions,
                        self.biomass_index,
                        self.biomass_function,
                    ) = read_mat_file(arg)
            except LookupError:
                print("An unexpected error occured when reading the input file.")
                sys.exit(1)

        elif isinstance(args, tuple):

            if len(arg) != 7:
                raise Exception(
                    "An unknown input format given to initialize a metabolic network object."
                )
            else:
                self.lb = arg[0]
                self.ub = arg[1]
                self.S = arg[2]
                self.metabolites = arg[3]
                self.reactons = arg[4]
                self.biomass_index = arg[5]
                self.biomass_function = arg[6]

                try:
                    if (
                        self.lb.size != self.ub.size
                        or self.lb.size != S.shape[1]
                        or self.metabolites.size != S.shape[0]
                        or self.reaction.size != S.shape[1]
                        or self.biomass_function.size != S.shape[1]
                        or (self.biomass_index < 0)
                        or (self.biomass_index > self.biomass_function.size)
                    ):
                        raise Exception(
                            "Wrong tuple format given to initialize a metabolic network object."
                        )
                except LookupError:
                    print(
                        "Wrong tuple format given to initialize a metabolic network object."
                    )
                    sys.exit(1)

        else:
            raise Exception(
                "An unknown input format given to initialize a metabolic network object."
            )

    def fva(self):
        """A member function to apply the FVA method on the metabolic network."""

        if self.parameters["fast_computations"]:
            (
                min_fluxes,
                max_fluxes,
                max_biomass_flux_vector,
                max_biomass_objective,
            ) = fast_fva(
                self.lb,
                self.ub,
                self.S,
                self.biomass_function,
                self.parameters["opt_percentage"],
            )
        else:
            (
                min_fluxes,
                max_fluxes,
                max_biomass_flux_vector,
                max_biomass_objective,
            ) = slow_fva(
                self.lb,
                self.ub,
                self.S,
                self.biomass_function,
                self.parameters["opt_percentage"],
            )

        return min_fluxes, max_fluxes, max_biomass_flux_vector, max_biomass_objective

    def fba(self):
        """A member function to apply the FBA method on the metabolic network."""

        if self.parameters["fast_computations"]:
            fba_res = fast_fba(self.lb, self.ub, self.S, self.biomass_function)
        else:
            fba_res = slow_fba(self.lb, self.ub, self.S, self.biomass_function)

        return fba_res

    @property
    def lb(self):
        return self.__lb

    @property
    def ub(self):
        return self.__ub

    @property
    def S(self):
        return self.__S

    @property
    def metabolites(self):
        return self.__metabolites

    @property
    def reactions(self):
        return self.__reactions

    @property
    def biomass_index(self):
        return self.__biomass_index

    @property
    def biomass_function(self):
        return self.__biomass_function

    @property
    def get_as_tuple(self):
        return (
            self.__lb,
            self.__ub,
            self.__S,
            self.__metabolites,
            self.__reactions,
            self.__biomass_index,
            self.__biomass_function,
        )

    def num_of_reactions(self):
        return len(self.reactions)

    @lb.setter
    def lb(self, value):
        self.__lb = value

    @ub.setter
    def ub(self, value):
        self.__ub = value

    @S.setter
    def S(self, value):
        self.__S = value

    @metabolites.setter
    def metabolites(self, value):
        self.__metabolites = value

    @reactions.setter
    def reactions(self, value):
        self.__reactions = value

    @biomass_index.setter
    def biomass_index(self, value):
        self.__biomass_index = value

    @biomass_function.setter
    def biomass_function(self, value):
        self.__biomass_function = value

    def set_fast_mode(self):

        try:
            import gurobipy

            self.parameters["fast_computations"] = True
        except ImportError as e:
            print("You have to install gurobi to use the fast computations.")
            self.parameters["fast_computations"] = False

    def set_slow_mode(self):

        self.parameters["fast_computations"] = False

    def set_nullspace_method(self, value):

        self.parameters["nullspace_method"] = value

    def set_opt_percentage(self, value):

        self.parameters["opt_percentage"] = value
