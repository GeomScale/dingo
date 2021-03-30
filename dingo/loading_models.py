# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import json
import scipy.io
import numpy as np


def read_json_file(input_file):
    """A Python function to Read a Bigg json file and returns,
    (a) lower/upper flux bounds
    (b) the stoichiometric matrix S (dense format)
    (c) the list of the metabolites
    (d) the list of reactions
    (e) the index of the biomass pseudoreaction
    (f) the objective function to maximize the biomass pseudoreaction

    Keyword arguments:
    input_file -- a json file that contains the information about a mettabolic network, for example see http://bigg.ucsd.edu/models
    """

    with open(input_file, "r") as f:

        distros_dict = json.load(f)

        reactions_list = distros_dict["reactions"]

        metabolites = []
        reactions = []

        for reaction in reactions_list:

            metabolites_dic = reaction["metabolites"]
            reaction_name = reaction["id"]
            reactions.append(reaction_name)

            for metabolite in metabolites_dic.keys():
                if metabolite not in metabolites:
                    metabolites.append(metabolite)

        list_of_reaction_lists = []
        vector_of_ubs = []
        vector_of_lbs = []

        for reaction in reactions_list:

            ub = float(reaction["upper_bound"])
            vector_of_ubs.append(ub)
            lb = float(reaction["lower_bound"])
            vector_of_lbs.append(lb)

            metabolites_dic = reaction["metabolites"]
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

    lb = np.asarray(lb, dtype="float")
    lb = np.ascontiguousarray(lb, dtype="float")

    ub = np.asarray(ub, dtype="float")
    ub = np.ascontiguousarray(ub, dtype="float")

    # The stoichiometric martrix S
    S = np.asarray(list_of_reaction_lists)
    S = np.transpose(S)

    S = np.asarray(S, dtype="float")
    S = np.ascontiguousarray(S, dtype="float")

    # Get biomass function if there
    biomass_function = np.zeros(S.shape[1])
    biomass_index = None
    for i in reactions:
        j = i.casefold()
        if "biom" in j:
            biomass_index = reactions.index(i)
            biomass_function[biomass_index] = 1

    return lb, ub, S, metabolites, reactions, biomass_index, biomass_function


def read_mat_file(input_file):
    """A Python function to Read a Bigg mat file and returns,
    (a) lower/upper flux bounds
    (b) the stoichiometric matrix S (dense format)
    (c) the list of the metabolites
    (d) the list of reactions
    (e) the index of the biomass pseudoreaction
    (f) the objective function to maximize the biomass pseudoreaction

    Keyword arguments:
    input_file -- a mat file that contains a MATLAB structure with the information about a mettabolic network, for example see http://bigg.ucsd.edu/models
    """

    data_from_mat = scipy.io.loadmat(input_file)

    species_name = ""
    for key in data_from_mat.keys():
        if key[0] != "_":
            species_name = key

    species = data_from_mat[species_name]
    list_of_lists = species.tolist()

    counter = 0

    for element in list_of_lists[0][0]:
        
        # position 0 corresponds to the Stoichiometric matrix
        if counter == 0:
            S = element
        # position 1 corresponds to the lower bounds
        if counter == 1:
            lb_tmp = element
            lb_tmp = lb_tmp.tolist()
        # position 2 corresponds to the upper bounds
        if counter == 2:
            ub_tmp = element
        # position 3 corresponds to the objective function
        if counter == 3:
            c_tmp = element
        # position 5 corresponds to the reactions
        if counter == 5:
            reactions_list = element.tolist()
            reactions = [reaction[0][0] for reaction in reactions_list]
        # position 6 corresponds to the metabolites
        if counter == 6:
            metabolites_list = element.tolist()
            metabolites = [metabolite[0][0] for metabolite in metabolites_list]
        
        counter += 1

    # Build function's output

    # lower and upper flux bounds
    ub = [i[0] for i in ub_tmp]
    lb = [x[0] for x in lb_tmp]
    biomass_function = [x[0] for x in c_tmp]

    lb = np.asarray(lb)
    ub = np.asarray(ub)
    biomass_function = np.asarray(biomass_function)

    biomass_function = np.asarray(biomass_function, dtype="float")
    biomass_function = np.ascontiguousarray(biomass_function, dtype="float")

    lb = np.asarray(lb, dtype="float")
    lb = np.ascontiguousarray(lb, dtype="float")

    ub = np.asarray(ub, dtype="float")
    ub = np.ascontiguousarray(ub, dtype="float")

    # The stoichiometric martrix S
    S = np.asarray(S, dtype="float")
    S = np.ascontiguousarray(S, dtype="float")

    # Get biomass index
    biomass_index = np.where(biomass_function == 1)
    biomass_index = biomass_index[0][0]

    return lb, ub, S, metabolites, reactions, biomass_index, biomass_function
