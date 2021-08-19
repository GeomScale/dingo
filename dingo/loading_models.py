# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import json, os
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


# This implementation works only for communities of 2 models 
# Once our method is validated, we will move on to implement it for more 
# by using the buildConqMatrix function
def getModelList(directory, format_type):

    """
    A Python function to get all the metabolic network files under a directory 
    and build a concatenated model and return:
    (a) lower/upper flux bounds
    (b) the stoichiometric matrix S (dense format)
    (c) the list of the metabolites
    (d) the list of reactions
    (e) the index of the biomass pseudoreaction
    (f) the objective function to maximize the biomass pseudoreaction

    Keyword arguments:
    directory -- directory where the metabolic network files of interest are located 
    """

    from dingo.MetabolicNetwork import MetabolicNetwork

    modelList = []
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            if format_type == "mat":
                model = MetabolicNetwork.from_mat(f)
            if format_type == "json":
                model = MetabolicNetwork.from_json(f)
            modelList.append(model)

    # At a a later point
    # Come back and replace the following section with a loop 
    # where conc_S will be like: 
    # S1   0   0   0 
    # 0    S2  0   0 
    # 0    0   S3  0 
    # 0    0   0   S4 ....   

    model_A = modelList[0]
    model_B = modelList[1]
    
    # Build concatenated stoichiometric matrix
    compl_1 = np.zeros((model_A.S.shape[0], model_B.S.shape[1]))
    compl_2 = np.zeros((model_B.S.shape[0], model_A.S.shape[1]))
    part_a  = np.concatenate((model_A.S, compl_1), axis=1)
    part_b  = np.concatenate((model_B.S, compl_2), axis=1)

    # Build concatenated biomass function 
    pair_biomass_function    = np.concatenate((model_A.S[:,model_A.biomass_index], model_B.S[:,model_B.biomass_index]), axis=0)
    pair_biomass_function    = pair_biomass_function.reshape((pair_biomass_function.shape[0], 1))
    conc_S  = np.concatenate((part_a, part_b), axis=0)
    conc_S  = np.concatenate((conc_S, pair_biomass_function), axis=1)
            
    # Get concatenated bounds
    conc_lb = np.concatenate((model_A.lb, model_B.lb), axis=0)
    conc_lb = np.append(conc_lb, 0.0)
    conc_ub = np.concatenate((model_A.ub, model_B.ub), axis=0)
    conc_ub = np.append(conc_ub, 1000.0)
    
    # Get concatenated reactions.. (including biomass overall)
    conc_reactions = model_A.reactions + model_B.reactions 
    conc_reactions.append("biomass_overall")
    conc_reactions = conc_reactions

    # .. and metabolites
    conc_metabolites = model_A.metabolites + model_B.metabolites

    # Overall biomass function info
    conc_biomass_function      = np.zeros((1,conc_S.shape[1] -1))
    conc_biomass_function = np.append(conc_biomass_function, 1.0)
    conc_biomass_index    = conc_S.shape[1] -1 

    return conc_lb, conc_ub, conc_S, conc_metabolites, conc_reactions, conc_biomass_index, conc_biomass_function, modelList