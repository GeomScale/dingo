# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import argparse


def dingo_args():
   parser = argparse.ArgumentParser()

   parser = argparse.ArgumentParser(
         description = 
         "dingo is a Python library for the analysis of \
         metabolic networks developed by the \
         GeomScale group - https://geomscale.github.io/ ",
         
         usage='%(prog)s [--help | -h] : help \n\n \
         1. by providing just your metabolic model: \n \
         python dingo.py -i my_model \n\n \
         2. or by asking for more: \n \
         python dingo.py -i my_model  -n 2000 -s gurobi \n \
         \n '
         )

   parser._action_groups.pop()

   required = parser.add_argument_group('required arguments')
   required.add_argument('--metabolic-network', '-i', 
      help = 'metabolic network as a .json or .mat file', 
      required=True,
      metavar = ''
      )

   optional = parser.add_argument_group('optional arguments')
   optional.add_argument('--sample-points', '-n', 
      help = 'number of sampling points to return', 
      required = False, 
      default = 1000,
      metavar = ''
      )

   optional.add_argument('--output-directory', '-o', 
      help = 'output directory for the dingo output', 
      required = False, 
      metavar = ''
      )

   optional.add_argument('--nullspace-method', '-m', 
      help = 'nullspace methods', 
      required = False, 
      default = 'sparsesvd',
      metavar = ''
      )

   optional.add_argument('--distribution', '-d', 
      help = 'choose among uniform, gaussian and exponential distribution for sampling the flux space of the metabolic network', 
      required = False, 
      default = 'uniform',
      metavar = ''
      )

   optional.add_argument('--solver', '-s', 
      help = 'if available you may set this parameter equal to "gurobi" to sample faster', 
      required = False, 
      default = 'linprog',
      metavar = ''
      )


   args = parser.parse_args()
   return args

