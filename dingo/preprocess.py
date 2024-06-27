
import cobra
from cobra.io import load_json_model
import cobra.manipulation
from collections import Counter


class PreProcess:
    
    def __init__(self, model):
        self.model = model
        
        
    def objective_function(model):
        
        objective = str(model.summary()._objective)
        objective = objective.split(" ")[1]
        return objective
    

    def zero_flux(model):
        
        tol = 1e-6
        
        fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0)
        zero_flux = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        zero_flux = zero_flux.index.tolist()
        
        return zero_flux
    
    
    def blocked(model):
        
        return cobra.flux_analysis.find_blocked_reactions(model)
    

    def metabolically_less_efficient(model):
        
        objective = PreProcess.objective_function(model)
    
        tol = 1e-6

        model.objective = objective
        fba_solution = model.optimize()

        wt_lower_bound = model.reactions.get_by_id(objective).lower_bound
        model.reactions.get_by_id(objective).lower_bound = fba_solution.objective_value

        fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.95)
        mle = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        mle = mle.index.tolist()
        
        model.reactions.get_by_id(objective).lower_bound = wt_lower_bound
        
        return mle
    
    
    def list_removed_reactions(model):

        remove_reactions = []

        blocked = PreProcess.blocked(model)
        mle = PreProcess.metabolically_less_efficient(model)
        zero_flux = PreProcess.zero_flux(model)

        remove_reactions = blocked+mle+zero_flux
        list_removed_reactions = list(set(remove_reactions))
        
        return list_removed_reactions
        

    def remove_model_reactions(model):
                
        removed_reactions_list = PreProcess.list_removed_reactions(model)
        
        for reaction in removed_reactions_list:
            model.reactions.get_by_id(reaction).lower_bound = 0
            model.reactions.get_by_id(reaction).upper_bound = 0
            
        return removed_reactions_list
            

    def possible_essential_reactions(model):
        
        tol = 1e-6
        removed_reactions_list = PreProcess.list_removed_reactions(model)  
    
        # find model reactions
        reactions_list = []
        
        for reaction in model.reactions:
            reaction_id = reaction.id
            reactions_list.append(reaction_id)
            
        remained_reactions = list((Counter(reactions_list)-Counter(removed_reactions_list)).elements())
   
        # find essential reactions
        essential_reactions_list = []
        essential_reactions = cobra.flux_analysis.find_essential_reactions(model)
        for reaction in essential_reactions:
            reaction_id = reaction.id
            essential_reactions_list.append(reaction_id)
            
        possible_essential = list((Counter(remained_reactions)-Counter(essential_reactions_list)).elements())
        print(possible_essential)
        
        
        PreProcess.remove_model_reactions(model)
        
        
        final_possible_essential = []

        for reaction in possible_essential:
            
            fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.9)
            disabled_before = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
            disabled_before = len(disabled_before.index.tolist())
            
            initial_lower = model.reactions.get_by_id(reaction).lower_bound
            initial_upper = model.reactions.get_by_id(reaction).upper_bound
                      
            model.reactions.get_by_id(reaction).lower_bound = 0
            model.reactions.get_by_id(reaction).upper_bound = 0
    
            fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.9)
            disabled_after = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
            disabled_after = len(disabled_after.index.tolist())
                                   
            if disabled_before != disabled_after:
                final_possible_essential.append(reaction)
                model.reactions.get_by_id(reaction).upper_bound = initial_upper
                model.reactions.get_by_id(reaction).lower_bound = initial_lower
            
        return final_possible_essential
        


model = load_json_model("../ext_data/e_coli_core.json")

fba_solution = model.optimize()
print(fba_solution.objective_value)

#possible_essentials = PreProcess.possible_essential_reactions(model)
#print(possible_essentials)

model.reactions.get_by_id("PFK").lower_bound = 0
model.reactions.get_by_id("PFK").upper_bound = 0


fba_solution = model.optimize()
print(fba_solution.objective_value)


