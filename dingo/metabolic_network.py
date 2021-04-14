# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import sys
from ding.loading_models import read_json_file, read_mat_file

class metabolic_network:

    def __init__(self, arg):

        self.opt_percentage = 100

        if isinstance(arg, str):

            if arg[-3:] != "mat" and arg[-4:] != "json":
                raise Exception("An unknown format file given.")
            try:
                if arg[-4:] == "json":
                    self.lb, self.ub, self.S, self.metabolites, self.reactons, self.biomass_index, self.biomass_function = read_json_file(arg)
                else:
                    self.lb, self.ub, self.S, self.metabolites, self.reactons, self.biomass_index, self.biomass_function = read_json_file(arg)
            except LookupError:
                print("An unexpected error occured when reading the input file.")
                sys.exit(1)
            
        elif isinstance(args[0], tuple):

            if (len(arg)!=7):
                raise Exception("An unknown input format given to initialize a metabolic network object.")
            else:
                self.lb = arg[0]
                self.ub = arg[1]
                self.S = arg[2]
                self.metabolites = arg[3]
                self.reactons = arg[4]
                self.biomass_index = arg[5]
                self.biomass_function = arg[6]

                try:
                    if (self.lb.size != self.ub.size or self.lb.size != S.shape[1] or self.metabolites.size != S.shape[0] or self.reaction.size != S.shape[1], or self.biomass_function.size != S.shape[1] or (self.biomass_index < 0 or self.biomass_index > self.biomass_function.size)):
                        raise Exception("Wrong tuple format given to initialize a metabolic network object.")
                except LookupError:
                    print("Wrong tuple format given to initialize a metabolic network object.")
                    sys.exit(1)

        else:
            raise Exception("An unknown input format given to initialize a metabolic network object.")
    
    def apply_fva(self):

        fva_res = fast_fva(self.lb, self.ub, self.S, self.biomass_function, self.opt_percentage)

        AA = fva_res[0] # for testing
        bb = fva_res[1] # for testing
        Aeqq = fva_res[2] # for testing
        beqq = fva_res[3] # for testing
        min_fluxes = fva_res[4]
        max_fluxes = fva_res[5]
        max_biomass_flux_vector = fva_res[6]
        max_biomass_objective = fva_res[7]

        return min_fluxes, max_fluxes, max_biomass_flux_vector, max_biomass_objective

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
        return self.__lb, self.__ub, self.__S, self.__metabolites, self.__reactions, self.__biomass_index, self.__biomass_function

    @data.setter
    def lb(self, value):
        self.__lb = value
    
    @data.setter
    def ub(self, value):
        self.__ub = value
    
    @data.setter
    def S(self, value):
        self.__S = value
    
    @data.setter
    def metabolites(self, value):
        self.__metabolites = value
    
    @data.setter
    def reactions(self, value):
        self.__reactions = value
    
    @data.setter
    def biomass_index(self, value):
        self.__biomass_index = value
    
    @data.setter
    def biomass_function(self, value):
        self.__biomass_function = value
