# Dingo : a python library for metabolic networks sampling and analysis
# Dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import json
import scipy.io
import numpy as np


## Read a Bigg file and get,
#     (a) lower/upper flux bounds
#     (b) the stoichiometric matrix S (dense format)
#     (c) the list of the metabolites
#     (d) the list of reactions

# The .json format case
def read_json_file(input_file):

   with open(input_file, 'r') as f:

      distros_dict = json.load(f)

      reactions_list = distros_dict['reactions']

      metabolites = []
      reactions = []

      for reaction in reactions_list:

         metabolites_dic = reaction['metabolites']
         reaction_name = reaction['id']
         reactions.append(reaction_name)

         for metabolite in metabolites_dic.keys():
            if metabolite not in metabolites:
               metabolites.append(metabolite)

      list_of_reaction_lists = []
      vector_of_ubs = []
      vector_of_lbs = []

      for reaction in reactions_list:

         ub = float(reaction['upper_bound']) ; vector_of_ubs.append(ub)
         lb = float(reaction['lower_bound']) ; vector_of_lbs.append(lb)

         metabolites_dic = reaction['metabolites']
         reaction_vector = []

         for term in range(len(metabolites)):

            metabolite = metabolites[term]

            if metabolite in metabolites_dic.keys():

               reaction_vector.append(metabolites_dic[metabolite])
            else:
               reaction_vector.append(0)

         list_of_reaction_lists.append(reaction_vector)

   # Build function's output; 
   
   # lower and upper flux bounds
   lb = np.asarray(vector_of_lbs)
   ub = np.asarray(vector_of_ubs)

   lb = np.asarray(lb, dtype = 'float')
   lb = np.ascontiguousarray(lb, dtype = 'float')

   ub = np.asarray(ub, dtype = 'float')
   ub = np.ascontiguousarray(ub, dtype = 'float')

   # The stoichiometric martrix S
   S = np.asarray(list_of_reaction_lists)
   S = np.transpose(S)

   S = np.asarray(S, dtype = 'float')
   S = np.ascontiguousarray(S, dtype = 'float')

   return lb, ub, S, metabolites, reactions


# The .mat format case
def read_mat_file(input_file):

   data_from_mat = scipy.io.loadmat(input_file)

   species_name = ''
   for key in data_from_mat.keys():
      if key[0] != "_":
         species_name = key

   species = data_from_mat[species_name]
   list_of_lists = species.tolist()

   counter = 0

   metabolites = []

   for element in list_of_lists[0][0]:

      if counter == 0:

         m =len(element)

         for i in element:

            metabolite = i[0][0]

            if metabolite not in metabolites:
               metabolites.append(metabolite)
      #position 7 corresponds to the reactions
      if counter == 7:
         reactions_list = element.tolist()
         reactions = [reaction[0][0] for reaction in reactions_list]
      #position 11 corresponds to the lower bounds
      if counter == 11:
         lb_tmp = element
         lb_tmp = lb_tmp.tolist()
      #position 12 corresponds to the upper bounds
      if counter == 12:
         ub_tmp = element
      #position 10 corresponds to the Stoichiometric matrix
      if counter == 10:
         Aeq = element

      counter += 1

   # Build function's output;

   # lower and upper flux bounds
   ub = [i[0] for i in ub_tmp]
   lb = [x[0] for x in lb_tmp]

   lb = np.asarray(lb)
   ub = np.asarray(ub)

   lb = np.asarray(lb, dtype = 'float')
   lb = np.ascontiguousarray(lb, dtype = 'float')

   ub = np.asarray(ub, dtype = 'float')
   ub = np.ascontiguousarray(ub, dtype = 'float')

   # The stoichiometric martrix S
   S = np.asarray(Aeq)

   return lb, ub, S, metabolites, reactions

