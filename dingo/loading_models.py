# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2021 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import json
import numpy as np
import cobra

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
    """A Python function based on the  to read a .mat file and returns,
    (a) lower/upper flux bounds
    (b) the stoichiometric matrix S (dense format)
    (c) the list of the metabolites
    (d) the list of reactions
    (e) the index of the biomass pseudoreaction
    (f) the objective function to maximize the biomass pseudoreaction

    Keyword arguments:
    input_file -- a mat file that contains a MATLAB structure with the information about a mettabolic network, for example see http://bigg.ucsd.edu/models
    """
    try: 
        cobra.io.load_matlab_model(  input_file )
    except:
        cobra_config = cobra.Configuration()
        cobra_config.solver = 'glpk'

    model = cobra.io.load_matlab_model( input_file )

    return (parse_cobra_model( model ))

def read_sbml_file(input_file):
    """A Python function, based on the cobra.io.read_sbml_model() function of cabrapy  
    and the extract_polytope() function of PolyRound 
    (https://gitlab.com/csb.ethz/PolyRound/-/blob/master/PolyRound/static_classes/parse_sbml_stoichiometry.py)
    to read an SBML file (.xml) and return:
    (a) lower/upper flux bounds
    (b) the stoichiometric matrix S (dense format)
    (c) the list of the metabolites
    (d) the list of reactions
    (e) the index of the biomass pseudoreaction
    (f) the objective function to maximize the biomass pseudoreaction

    Keyword arguments:
    input_file -- a xml file that contains an SBML  model with the information about a mettabolic network, for example see: 
    https://github.com/VirtualMetabolicHuman/AGORA/blob/master/CurrentVersion/AGORA_1_03/AGORA_1_03_sbml/Abiotrophia_defectiva_ATCC_49176.xml
    """
    try: 
        cobra.io.read_sbml_model(  input_file )
    except:
        cobra_config = cobra.Configuration()
        cobra_config.solver = 'glpk'

    model = cobra.io.read_sbml_model( input_file )

    return (parse_cobra_model( model ))

def parse_cobra_model(cobra_model):

    inf_bound=1e5

    metabolites = [ metabolite.id for metabolite in cobra_model.metabolites ]
    reactions = [ reaction.id for reaction in cobra_model.reactions ]

    S = cobra.util.array.create_stoichiometric_matrix(cobra_model)

    lb  = []
    ub = []
    biomass_function = np.zeros( len(cobra_model.reactions) )

    for index, reaction in enumerate(cobra_model.reactions):

        if reaction.objective_coefficient==1:
            biomass_index = index
            biomass_function[index] = 1

        if reaction.bounds[0] == float("-inf"):
            lb.append( -inf_bound )
        else:
            lb.append( reaction.bounds[0] )

        if reaction.bounds[1] == float("inf"):
            ub.append( inf_bound )
        else:
            ub.append( reaction.bounds[1] )

    lb = np.asarray(lb)
    ub = np.asarray(ub)

    biomass_function = np.asarray(biomass_function)
    biomass_function = np.asarray(biomass_function, dtype="float")
    biomass_function = np.ascontiguousarray(biomass_function, dtype="float")

    lb = np.asarray(lb, dtype="float")
    lb = np.ascontiguousarray(lb, dtype="float")

    ub = np.asarray(ub, dtype="float")
    ub = np.ascontiguousarray(ub, dtype="float")

    return lb, ub, S, metabolites, reactions, biomass_index, biomass_function
