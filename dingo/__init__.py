# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import sys
import os
import pickle
from dingo.loading_models import read_json_file
from dingo.nullspace import nullspace_dense, nullspace_sparse
from dingo.scaling import gmscale
from dingo.utils import (
    apply_scaling,
    remove_almost_redundant_facets,
    map_samples_to_steady_states,
    get_matrices_of_low_dim_polytope,
    get_matrices_of_full_dim_polytope,
)
from dingo.illustrations import (
    plot_copula,
    plot_histogram,
)
from dingo.parser import dingo_args
from dingo.MetabolicNetwork import MetabolicNetwork
from dingo.PolytopeSampler import PolytopeSampler

from dingo.pyoptinterface_based_impl import fba, fva, inner_ball, remove_redundant_facets, set_default_solver

from volestipy import HPolytope


def get_name(args_network):

    position = [pos for pos, char in enumerate(args_network) if char == "/"]

    if args_network[-4:] == "json":
        if position == []:
            name = args_network[0:-5]
        else:
            name = args_network[(position[-1] + 1) : -5]
    elif args_network[-3:] == "mat":
        if position == []:
            name = args_network[0:-4]
        else:
            name = args_network[(position[-1] + 1) : -4]

    return name


def dingo_main():
    """A function that (a) reads the inputs using argparse package, (b) calls the proper dingo pipeline
    and (c) saves the outputs using pickle package
    """

    args = dingo_args()

    if args.metabolic_network is None and args.polytope is None and not args.histogram:
        raise Exception(
            "You have to give as input either a model or a polytope derived from a model."
        )

    if args.metabolic_network is None and ((args.fva) or (args.fba)):
        raise Exception("You have to give as input a model to apply FVA or FBA method.")

    if args.output_directory == None:
        output_path_dir = os.getcwd()
    else:
        output_path_dir = args.output_directory

    if os.path.isdir(output_path_dir) == False:
        os.mkdir(output_path_dir)

    # Move to the output directory
    os.chdir(output_path_dir)
    
    set_default_solver(args.solver)

    if args.model_name is None:
        if args.metabolic_network is not None:
            name = get_name(args.metabolic_network)
    else:
        name = args.model_name

    if args.histogram:

        if args.steady_states is None:
            raise Exception(
                "A path to a pickle file that contains steady states of the model has to be given."
            )

        if args.metabolites_reactions is None:
            raise Exception(
                "A path to a pickle file that contains the names of the metabolites and the reactions of the model has to be given."
            )

        if int(args.reaction_index) <= 0:
            raise Exception("The index of the reaction has to be a positive integer.")

        file = open(args.steady_states, "rb")
        steady_states = pickle.load(file)
        file.close()

        file = open(args.metabolites_reactions, "rb")
        model = pickle.load(file)
        file.close()

        reactions = model.reactions

        if int(args.reaction_index) > len(reactions):
            raise Exception(
                "The index of the reaction has not to be exceed the number of reactions."
            )

        if int(args.n_bins) <= 0:
            raise Exception("The number of bins has to be a positive integer.")

        plot_histogram(
            steady_states[int(args.reaction_index) - 1],
            reactions[int(args.reaction_index) - 1],
            int(args.n_bins),
        )

    elif args.fva:

        if args.metabolic_network[-4:] == "json":
            model = MetabolicNetwork.fom_json(args.metabolic_network)
        elif args.metabolic_network[-3:] == "mat":
            model = MetabolicNetwork.fom_mat(args.metabolic_network)
        else:
            raise Exception("An unknown format file given.")

        model.set_solver(args.solver)
        
        result_obj = model.fva()

        with open("dingo_fva_" + name + ".pckl", "wb") as dingo_fva_file:
            pickle.dump(result_obj, dingo_fva_file)

    elif args.fba:

        if args.metabolic_network[-4:] == "json":
            model = MetabolicNetwork.fom_json(args.metabolic_network)
        elif args.metabolic_network[-3:] == "mat":
            model = MetabolicNetwork.fom_mat(args.metabolic_network)
        else:
            raise Exception("An unknown format file given.")
        
        model.set_solver(args.solver)

        result_obj = model.fba()

        with open("dingo_fba_" + name + ".pckl", "wb") as dingo_fba_file:
            pickle.dump(result_obj, dingo_fba_file)

    elif args.metabolic_network is not None:

        if args.metabolic_network[-4:] == "json":
            model = MetabolicNetwork.fom_json(args.metabolic_network)
        elif args.metabolic_network[-3:] == "mat":
            model = MetabolicNetwork.fom_mat(args.metabolic_network)
        else:
            raise Exception("An unknown format file given.")

        sampler = PolytopeSampler(model)

        if args.preprocess_only:

            sampler.get_polytope()

            polytope_info = (
                sampler,
                name,
            )

            with open("dingo_model_" + name + ".pckl", "wb") as dingo_model_file:
                pickle.dump(model, dingo_model_file)

            with open(
                "dingo_polytope_sampler_" + name + ".pckl", "wb"
            ) as dingo_polytope_file:
                pickle.dump(polytope_info, dingo_polytope_file)

        else:

            steady_states = sampler.generate_steady_states(
                int(args.effective_sample_size),
                args.psrf_check,
                args.parallel_mmcs,
                int(args.num_threads),
            )

            polytope_info = (
                sampler,
                name,
            )

            with open("dingo_model_" + name + ".pckl", "wb") as dingo_model_file:
                pickle.dump(model, dingo_model_file)

            with open(
                "dingo_polytope_sampler_" + name + ".pckl", "wb"
            ) as dingo_polytope_file:
                pickle.dump(polytope_info, dingo_polytope_file)

            with open(
                "dingo_steady_states_" + name + ".pckl", "wb"
            ) as dingo_steadystates_file:
                pickle.dump(steady_states, dingo_steadystates_file)

    else:

        file = open(args.polytope, "rb")
        input_obj = pickle.load(file)
        file.close()
        sampler = input_obj[0]

        if isinstance(sampler, PolytopeSampler):

            steady_states = sampler.generate_steady_states(
                int(args.effective_sample_size),
                args.psrf_check,
                args.parallel_mmcs,
                int(args.num_threads),
            )

        else:
            raise Exception("The input file has to be generated by dingo package.")

        if args.model_name is None:
            name = input_obj[-1]

        polytope_info = (
            sampler,
            name,
        )

        with open(
            "dingo_polytope_sampler" + name + "_improved.pckl", "wb"
        ) as dingo_polytope_file:
            pickle.dump(polytope_info, dingo_polytope_file)

        with open("dingo_steady_states_" + name + ".pckl", "wb") as dingo_network_file:
            pickle.dump(steady_states, dingo_network_file)


if __name__ == "__main__":

    dingo_main()
