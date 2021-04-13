# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import sys
import os
import pickle
from dingo.fva import slow_fva
from dingo.fba import slow_fba
from dingo.loading_models import read_json_file
from dingo.inner_ball import slow_inner_ball
from dingo.nullspace import nullspace_dense, nullspace_sparse
from dingo.scaling import (
    gmscale,
    apply_scaling,
    remove_almost_redundant_facets,
    map_samples_to_steady_states,
)
from dingo.parser import dingo_args
from dingo.pipelines import (
    from_model_to_steady_states_pipeline,
    from_polytope_to_steady_states_pipeline,
    fva_pipeline,
    fba_pipeline,
    plot_histogram,
)

try:
    import gurobipy
    from dingo.gurobi_based_implementations import fast_fba, fast_fva, fast_inner_ball
except ImportError as e:
    pass

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
    print(name)
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

        steady_states = steady_states[0]

        file = open(args.metabolites_reactions, "rb")
        meta_reactions = pickle.load(file)
        file.close()

        reactions = meta_reactions[1]

        plot_histogram(
            steady_states[int(args.reaction_index) - 1],
            reactions[int(args.reaction_index) - 1],
            int(args.n_bins),
        )

    elif args.fva:

        result_obj = fva_pipeline(args)
        result_obj = result_obj[4:]

        with open("dingo_fva_" + name, "wb") as dingo_fva_file:
            pickle.dump(result_obj, dingo_fva_file)

    elif args.fba:

        result_obj = fba_pipeline(args)

        with open("dingo_fba_" + name, "wb") as dingo_fba_file:
            pickle.dump(result_obj, dingo_fba_file)

    elif args.metabolic_network is not None:

        result_obj = from_model_to_steady_states_pipeline(args)

        if args.preprocess_only:

            polytope_info = result_obj[:4]
            network_info = result_obj[4:6]
            metabol_reaction = result_obj[6:]

            polytope_info = polytope_info + (name,)

            with open("dingo_polytope_" + name, "wb") as dingo_polytope_file:
                pickle.dump(polytope_info, dingo_polytope_file)

            with open(
                "dingo_metabolites_reactions_" + name, "wb"
            ) as dingo_metabolreactions_file:
                pickle.dump(metabol_reaction, dingo_metabolreactions_file)

            with open("dingo_minmax_fluxes_" + name, "wb") as dingo_network_file:
                pickle.dump(network_info, dingo_network_file)

        else:

            polytope_info = result_obj[:7]
            network_info = result_obj[7:9]
            metabol_reaction = result_obj[9:11]
            steadystates = result_obj[11:]

            polytope_info = polytope_info + (name,)

            with open("dingo_rounded_polytope_" + name, "wb") as dingo_polytope_file:
                pickle.dump(polytope_info, dingo_polytope_file)

            with open("dingo_minmax_fluxes_" + name, "wb") as dingo_network_file:
                pickle.dump(network_info, dingo_network_file)

            with open(
                "dingo_metabolites_reactions_" + name, "wb"
            ) as dingo_metabolreactions_file:
                pickle.dump(metabol_reaction, dingo_metabolreactions_file)

            with open("dingo_steady_states_" + name, "wb") as dingo_steadystates_file:
                pickle.dump(steadystates, dingo_steadystates_file)

    else:

        file = open(args.polytope, "rb")
        object_file = pickle.load(file)
        file.close()

        if len(object_file) == 5:
            result_obj = from_polytope_to_steady_states_pipeline(
                args, object_file[0], object_file[1], object_file[2], object_file[3]
            )
        elif len(object_file) == 8:
            result_obj = from_polytope_to_steady_states_pipeline(
                args,
                object_file[0],
                object_file[1],
                object_file[2],
                object_file[3],
                object_file[4],
                object_file[5],
            )
        else:
            raise Exception("The input file has to be generated by dingo package.")

        polytope_info = result_obj[:7]
        network_info = result_obj[7:]

        if args.model_name is None:
            name = object_file[-1]

        polytope_info = polytope_info + (name,)

        with open("dingo_double_rounded_polytope_" + name, "wb") as dingo_polytope_file:
            pickle.dump(polytope_info, dingo_polytope_file)

        with open("dingo_steady_states_" + name, "wb") as dingo_network_file:
            pickle.dump(network_info, dingo_network_file)


if __name__ == "__main__":

    dingo_main()
