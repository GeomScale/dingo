
import cobra
from cobra.io import load_json_model
import cobra.manipulation


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
        
        list_removed_reactions = PreProcess.list_removed_reactions(model)
        
        for reaction in list_removed_reactions:
            model.reactions.get_by_id(reaction).lower_bound = 0
            model.reactions.get_by_id(reaction).upper_bound = 0


model = load_json_model("../ext_data/e_coli_core.json")

#fba_solution = model.optimize()
#print(fba_solution.objective_value)

new_model = PreProcess.remove_model_reactions(model)

#fba_solution = model.optimize()
#print(fba_solution.objective_value)
