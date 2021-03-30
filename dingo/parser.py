# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import argparse


def dingo_args():
    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser(
        description="dingo is a Python library for the analysis of \
         metabolic networks developed by the \
         GeomScale group - https://geomscale.github.io/ ",
        usage="%(prog)s [--help | -h] : help \n\n \
         1. by providing just your metabolic model: \n \
         python dingo.py -i my_model \n\n \
         2. or by asking for more: \n \
         python dingo.py -i my_model  -n 2000 -s gurobi \n \
         \n ",
    )

    parser._action_groups.pop()

    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "--metabolic_network",
        "-i",
        help="the path to a metabolic network as a .json or a .mat file.",
        required=True,
        metavar="",
    )

    optional = parser.add_argument_group("optional arguments")
    optional.add_argument(
        "--effective_sample_size",
        "-n",
        help="the minimum effective sample size per marginal of the sample that the Multiphase Monte Carlo Sampling algorithm will return. The default value is 1000.",
        required=False,
        default=1000,
        metavar="",
    )

    optional.add_argument(
        "--output_directory",
        "-o",
        help="the output directory for the dingo output",
        required=False,
        metavar="",
    )

    optional.add_argument(
        "--nullspace",
        "-null",
        help="the method to compute the right nullspace of the stoichiometric matrix. Choose between `sparseQR` and `scipy`. The default method is `sparseQR`.",
        required=False,
        default="sparseQR",
        metavar="",
    )

    optional.add_argument(
        "--psrf_check",
        "-psrf",
        help="a boolean flag to request psrf < 1.1 for each marginal of the sample that the Multiphase Monte Carlo Sampling algorithm will return. The default value is `False`.",
        required=False,
        default=False,
        metavar="",
    )

    optional.add_argument(
        "--parallel_mmcs",
        "-pmmcs",
        help="a boolean flag to request sampling with parallel Multiphase Monte Carlo Sampling algorithm. The default value is `false`.",
        required=False,
        default=False,
        metavar="",
    )

    optional.add_argument(
        "--num_threads",
        "-nt",
        help="the number of threads to be used in parallel Multiphase Monte Carlo Sampling algorithm. The default number is 2.",
        required=False,
        default=2,
        metavar="",
    )

    optional.add_argument(
        "--distribution",
        "-d",
        help="the distribution to sample from the flux space of the metabolic network. Choose among `uniform`, `gaussian` and `exponential` distribution. The default value is `uniform`.",
        required=False,
        default="uniform",
        metavar="",
    )

    optional.add_argument(
        "--solver",
        "-s",
        help="the solver to use for the linear programs. Choose between `scipy` and `gurobi` (faster computations --- it needs a licence). The default value is `scipy`.",
        required=False,
        default="scipy",
        metavar="",
    )

    args = parser.parse_args()
    return args
