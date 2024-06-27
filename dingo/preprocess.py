
import cobra
from cobra.io import load_json_model
import cobra.manipulation
from collections import Counter


class PreProcess:
    
    def __init__(self, model):
        self.model = model
        self.objective = self.objective_function()
        self.zero_flux_reactions = self.zero_flux()
        self.blocked_reactions = self.blocked()
        self.mle_reactions = self.metabolically_less_efficient()
        self.removed_reactions = self.list_removed_reactions()
        self.essential_reactions = self.possible_essential_reactions()
        
        
    def objective_function(self):
        
        objective = str(self.model.summary()._objective)
        objective = objective.split(" ")[1]
        
        self.objective = objective
        return self.objective
    

    def zero_flux(self):
        
        tol = 1e-6
        
        fva = cobra.flux_analysis.flux_variability_analysis(self.model, fraction_of_optimum=0.9)
        zero_flux = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        zero_flux_reactions = zero_flux.index.tolist()
        
        self.zero_flux_reactions = zero_flux_reactions
        return self.zero_flux_reactions
    
    
    def blocked(self):
        
        blocked_reactions = cobra.flux_analysis.find_blocked_reactions(self.model)
        self.blocked_reactions =  blocked_reactions
        return self.blocked_reactions


    # to add documentation comments
    def metabolically_less_efficient(self):
            
        tol = 1e-6

        self.model.objective = self.objective
        fba_solution = self.model.optimize()

        wt_lower_bound = self.model.reactions.get_by_id(self.objective).lower_bound
        self.model.reactions.get_by_id(self.objective).lower_bound = fba_solution.objective_value

        fva = cobra.flux_analysis.flux_variability_analysis(self.model, fraction_of_optimum=0.95)
        mle = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        mle = mle.index.tolist()
        
        self.model.reactions.get_by_id(self.objective).lower_bound = wt_lower_bound
        
        self.mle_reactions = mle
        return self.mle_reactions
    
    
    def list_removed_reactions(self):

        remove_reactions = []

        remove_reactions = self.blocked_reactions + self.mle_reactions + self.zero_flux_reactions
        list_removed_reactions = list(set(remove_reactions))
        
        self.removed_reactions = list_removed_reactions
        return self.removed_reactions
        

    def remove_model_reactions(self):
                        
        for reaction in self.removed_reactions:
            self.model.reactions.get_by_id(reaction).lower_bound = 0
            self.model.reactions.get_by_id(reaction).upper_bound = 0
            
        return self.model
            

    def possible_essential_reactions(self):
            
        # find model reactions
        reactions_list = []
        
        for reaction in self.model.reactions:
            reaction_id = reaction.id
            reactions_list.append(reaction_id)
            
        remained_reactions = list((Counter(reactions_list)-Counter(self.removed_reactions)).elements())
   
        # find essential reactions
        essential_reactions_list = []
        essential_reactions = cobra.flux_analysis.find_essential_reactions(self.model)
        for reaction in essential_reactions:
            reaction_id = reaction.id
            essential_reactions_list.append(reaction_id)
            
        possible_essential = list((Counter(remained_reactions)-Counter(essential_reactions_list)).elements())
        
        
        self.remove_model_reactions()
        
        
        for reaction in possible_essential:
                        
            initial_lower = self.model.reactions.get_by_id(reaction).lower_bound
            initial_upper = self.model.reactions.get_by_id(reaction).upper_bound
            
            fba_solution_before = self.model.optimize().objective_value

            self.model.reactions.get_by_id(reaction).lower_bound = 0.01 * initial_lower
            self.model.reactions.get_by_id(reaction).upper_bound = 0.01 * initial_upper
    
            fba_solution_after = self.model.optimize().objective_value
                                   
            if fba_solution_after < fba_solution_before:
                if abs(fba_solution_before)-abs(fba_solution_after) > 0.3:
                    essential_reactions_list.append(reaction)
                    
                    self.model.reactions.get_by_id(reaction).upper_bound = initial_upper
                    self.model.reactions.get_by_id(reaction).lower_bound = initial_lower
            
            
        self.essential_reactions = essential_reactions_list
        return self.essential_reactions
        
        


model = load_json_model("../ext_data/e_coli_core.json")

fba_solution = model.optimize()
print(fba_solution.objective_value)

obj = PreProcess(model)

print(len(obj.essential_reactions))

fba_solution = model.optimize()
print(fba_solution.objective_value)


