
from dingo import MetabolicNetwork
import numpy as np
import pandas as pd
import cobra
from cobra.io import load_json_model


# function to load and benefit from both dingo's and cobrapy's model attributes
def load_model(model):
    dingo_model = MetabolicNetwork.from_json(model)
    cobra_model = load_json_model(model)
    return dingo_model, cobra_model

dingo_model = load_model('../ext_data/e_coli_core.json')[0]
cobra_model = load_model('../ext_data/e_coli_core.json')[1]


# function to get a list of all reaction ids
def reactions_ids(dingo_model):
    reactions = dingo_model.reactions
    return reactions

reactions = reactions_ids(dingo_model)


# function to find all unidirected reactions (-->)
def unidirected_reactions(cobra_model):
    uni_reactions = []
    for reaction in reactions:
        try:
            str = cobra_model.reactions.get_by_id(reaction).reaction
        except:
            print("failed to find" , reaction)
        
        # continue if you find this arrow and not this <=>
        try:
            str.index("-->")
            uni_reactions.append(reaction)   
        except:
            pass
        
    return uni_reactions

uni_reactions = unidirected_reactions(cobra_model)


# function to find metabolites participating in possible lumped reactions
def lumped_metabolites(dingo_model):

    stoichiometric_matrix = dingo_model.S
    metabolites = dingo_model.metabolites
    
    positives = stoichiometric_matrix > 0
    negatives = stoichiometric_matrix < 0
    
    exactly_one_positive = np.sum(positives, axis=1) == 1
    exactly_one_negative = np.sum(negatives, axis=1) == 1
    
    # get boolean type array of rows where only one positive and one negative coefficient exists
    matching_rows = np.logical_and(exactly_one_positive, exactly_one_negative)
    
    # get the corresponding rows indeces
    lumped_met_index = np.where(matching_rows)[0]
    
    # get the corresponding metabolites ids
    lumped_metabolites = []
    for met in lumped_met_index:
        lumped_metabolites.append(metabolites[met])
        
    # get the corresponding arrays rows
    lumped_met_rows = stoichiometric_matrix[lumped_met_index]
    
    return lumped_metabolites, lumped_met_rows

lumped_metabolites = lumped_metabolites(dingo_model)


# function to find possible lumped reactions based on lumped metabolites
def lumped_reactions(dingo_model):
    
    reactions = dingo_model.reactions
    
    # find reactions where lumped metabolites are substrates
    columns_with_negatives = np.any(lumped_metabolites[1] < 0, axis=0)
    lumped_reactions = [reactions[i] for i, has_negative in enumerate(columns_with_negatives) if has_negative]
    
    return lumped_reactions

lumped_reactions = lumped_reactions(dingo_model)


# function to match lumped reaction-substrate-product
def lumped_reactions_substrate_product(cobra_model, lumped_reactions):
        
    # define some cofactors- this list must be updated
    cofactors = [
        "coa_c",
        "atp_c",
        "amp_c",
        "adp_c",
        "nad_c",
        "nadh_c",
        "nadp_c",
        "nadph_c",
        "h2o_c",
        "h_c",
        "accoa_c",
        "h_e",
        "co2_c",
        "o2_c"]
    
    substrates = []
    products = []
    final_lumped_reactions = []

    for lumped in lumped_reactions:
        str = cobra_model.reactions.get_by_id(lumped).reaction 
        try:
            str.index("-->")
            split = str.split(" ")
            index = split.index("-->")
            for substrate in split[0:index]:
                # approach to remove coeeficients (1.0, 2.0)
                try:
                    float(substrate)
                except:
                    if substrate not in cofactors and (substrate != "+"):   
                        substrates.append(substrate)
            
            for product in split[index+1:]:
                try:
                    float(product)
                except:   
                    if (product not in cofactors) and (product != "+"):
                        products.append(product)
                        
            final_lumped_reactions.append(lumped)
        except:
            pass
        
    # find which product becomes substrate in another reaction
    for i in range(len(products)):
        if products[i] in substrates:
            j = substrates.index(products[i])
            print(final_lumped_reactions[i], "can be lumped with", final_lumped_reactions[j])
                    
reactions_substrates_products = lumped_reactions_substrate_product(cobra_model, lumped_reactions)