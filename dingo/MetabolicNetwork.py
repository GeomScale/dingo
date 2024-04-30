# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Vissarion Fisikopoulos

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import pandas as pd
import sys
from typing import Dict
import cobra
from dingo.loading_models import read_json_file, read_mat_file, read_sbml_file, parse_cobra_model
from dingo.fva import slow_fva
from dingo.fba import slow_fba
import logging
logger = logging.getLogger(__name__)

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
            self._parameters["tol"] = 1e-06
        except ImportError as e:
            self._parameters["fast_computations"] = False
            self._parameters["tol"] = 1e-03

        if len(tuple_args) != 12:
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
        self._medium = tuple_args[7]
        self._medium_indices = tuple_args[8]
        self._exchanges = tuple_args[9]
        self._reactions_map = tuple_args[10]
        self._metabolites_map = tuple_args[11]

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
        if (not isinstance(arg, str)) and ((arg[-3:] == "xml") or (arg[-4] == "sbml")):
            raise Exception(
                "An unknown input format given to initialize a metabolic network object."
            )

        return cls(read_sbml_file(arg))

    @classmethod
    def from_cobra_model(cls, arg):
        if (not isinstance(arg, cobra.core.model.Model)):
            raise Exception(
                "An unknown input format given to initialize a metabolic network object."
            )

        return cls(parse_cobra_model(arg))

    def _fva(self):
        """A member function to apply the FVA method on the metabolic network."""

        if self._parameters["fast_computations"]:
            min_fluxes, max_fluxes, max_biomass_flux_vector, max_biomass_objective = fast_fva(
                self._lb,
                self._ub,
                self._S,
                self._biomass_function,
                self._parameters["opt_percentage"],
            )
        else:
            min_fluxes, max_fluxes, max_biomass_flux_vector, max_biomass_objective = slow_fva(
                self._lb,
                self._ub,
                self._S,
                self._biomass_function,
                self._parameters["opt_percentage"],
            )
        self._min_fluxes = min_fluxes
        self._max_fluxes = max_fluxes
        self._opt_vector = max_biomass_flux_vector
        self._opt_value = max_biomass_objective
        return min_fluxes, max_fluxes, max_biomass_flux_vector, max_biomass_objective

    def _fba(self):
        """A member function to apply the FBA method on the metabolic network."""

        if self._parameters["fast_computations"]:
            opt_vector, opt_value = fast_fba(self._lb, self._ub, self._S, self._biomass_function)
        else:
            opt_vector, opt_value = slow_fba(self._lb, self._ub, self._S, self._biomass_function)
        self._opt_vector = opt_vector
        self._opt_value = opt_value
        return opt_vector, opt_value

    def fba(self):
        if not hasattr(self, '_opt_vector'):
            self._fba()
        fba_df = pd.DataFrame({'fluxes': self._opt_vector}, index=self._reactions)
        return fba_df

    def fva(self):
        if not hasattr(self, '_min_fluxes'):
            self._fva()
        fva_df = pd.DataFrame({'minimum': self._min_fluxes, 'maximum': self._max_fluxes}, index=self._reactions)
        return fva_df

    # Descriptors
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
    def medium(self):
        return self._medium

    @property
    def exchanges(self):
        return self._exchanges

    @property
    def parameters(self):
        return self._parameters

    @property
    def reactions_map(self):
        return self._reactions_map

    @property
    def metabolites_map(self):
        return self._metabolites_map

    @property
    def opt_value(self):
        return self._opt_value

    @property
    def opt_vector(self):
        return self._opt_vector

    @property
    def min_fluxes(self):
        return self._min_fluxes

    @property
    def max_fluxes(self):
        return self._max_fluxes

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
            self._medium,
            self._exchanges,
            self._reactions_map,
            self._metabolites_map
        )

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

    @medium.setter
    def medium(self, medium: Dict[str, float]) -> None:
        """Set the constraints on the model exchanges.

        `model.medium` returns a dictionary of the bounds for each of the
        boundary reactions, in the form of `{rxn_id: rxn_bound}`, where `rxn_bound`
        specifies the absolute value of the bound in direction of metabolite
        creation (i.e., lower_bound for `met <--`, upper_bound for `met -->`)

        Parameters
        ----------
        medium: dict
            The medium to initialize. medium should be a dictionary defining
            `{rxn_id: bound}` pairs.
        """

        def set_active_bound(reaction: str, reac_index: int, bound: float) -> None:
            """Set active bound.

            Parameters
            ----------
            reaction: cobra.Reaction
                Reaction to set
            bound: float
                Value to set bound to. The bound is reversed and set as lower bound
                if reaction has reactants (metabolites that are consumed). If reaction
                has reactants, it seems the upper bound won't be set.
            """
            if any(x < 0 for x in  list(self._S[:, reac_index])):
                self._lb[reac_index] = -bound
            elif any(x > 0 for x in  list(self._S[:, reac_index])):
                self._ub[reac_index] = bound

        # Set the given media bounds
        media_rxns = []
        exchange_rxns = frozenset(self.exchanges)
        for rxn_id, rxn_bound in medium.items():
            if rxn_id not in exchange_rxns:
                logger.warning(
                    f"{rxn_id} does not seem to be an an exchange reaction. "
                    f"Applying bounds anyway."
                )
            media_rxns.append(rxn_id)

            reac_index = self._reactions.index(rxn_id)

            set_active_bound(rxn_id, reac_index, rxn_bound)

        frozen_media_rxns = frozenset(media_rxns)

        # Turn off reactions not present in media
        for rxn_id in exchange_rxns - frozen_media_rxns:
            """
            is_export for us, needs to check on the S
            order reactions to their lb and ub
            """
            # is_export = rxn.reactants and not rxn.products
            reac_index = self._reactions.index(rxn_id)
            products = np.any(self._S[:,reac_index] > 0)
            reactants_exist = np.any(self._S[:,reac_index] < 0)
            is_export = True if not products and reactants_exist else False
            set_active_bound(
                rxn_id, reac_index, min(0.0, -self._lb[reac_index] if is_export else self._ub[reac_index])
            )
        self._medium = medium
        self._opt_vector = None

    @reactions_map.setter
    def reactions_map(self, value):
        self._reactions_map = value

    @metabolites_map.setter
    def metabolites_map(self, value):
        self._metabolites_map = value

    @min_fluxes.setter
    def min_fluxes(self, value):
        self._min_fluxes = value

    @max_fluxes.setter
    def max_fluxes(self, value):
        self._max_fluxes = value

    @opt_value.setter
    def opt_value(self, value):
        self._opt_value = value

    @opt_vector.setter
    def opt_vector(self, value):
        self._opt_vector = value


    # Attributes
    def num_of_reactions(self):
        return len(self._reactions)

    def num_of_metabolites(self):
        return len(self._metabolites)

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
