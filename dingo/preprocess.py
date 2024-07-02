
import cobra
from cobra.io import load_json_model
import cobra.manipulation
from collections import Counter


class PreProcess:
    
    def __init__(self, model):
        self.model = model
        self.objective = self.objective_function()
        self.initial_reactions = self.initial()
        self.essential_reactions = self.essentials()
        self.zero_flux_reactions = self.zero_flux()
        self.blocked_reactions = self.blocked()
        self.mle_reactions = self.metabolically_less_efficient()
        self.removed_reactions = []

        
    def objective_function(self):
        """
        A function  used to find the objective function of a model
        """
        
        objective = str(self.model.summary()._objective)
        objective = objective.split(" ")[1]
        
        self.objective = objective
        return self.objective


    def initial(self):
        
        self.initial_reactions = []
        
        for reaction in self.model.reactions:
            reaction_id = reaction.id
            self.initial_reactions.append(reaction_id)
            
        return self.initial_reactions
    
    
    def reaction_bounds_dictionary(self):
        
        self.reaction_bounds_dict = {
                  'reaction': (0, 100)
                  }
                
        for reaction_id in self.initial_reactions:
            bounds = self.model.reactions.get_by_id(reaction_id).bounds
            self.reaction_bounds_dict[reaction_id] = bounds


    def essentials(self):
        
        self.essential_reactions = []
        essential_reactions = cobra.flux_analysis.find_essential_reactions(self.model)
        for reaction in essential_reactions:
            reaction_id = reaction.id
            self.essential_reactions.append(reaction_id)
            
        return self.essential_reactions
    

    def zero_flux(self):
        """
        A function used to find zero-flux reactions.
        These reactions are the ones that have a flux equaled to 0
        when running FVA analysis with fraction of optimum set to 90%
        """
        
        tol = 1e-6
        
        fva = cobra.flux_analysis.flux_variability_analysis(self.model, fraction_of_optimum=0.9)
        zero_flux = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        zero_flux_reactions = zero_flux.index.tolist()
        
        self.zero_flux_reactions = zero_flux_reactions
        return self.zero_flux_reactions
    
    
    def blocked(self):
        """
        A function used to find blocked reactions.
        These reactions can not have any flux other than 0
        """
        
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
    
    
    def remove_model_reactions(self):
                                
        for reaction in self.removed_reactions:
            self.model.reactions.get_by_id(reaction).lower_bound = 0
            self.model.reactions.get_by_id(reaction).upper_bound = 0
            
        return self.model
        
            
    def removed(self, extend):
        
        tol = 1e-6
        
        if extend != 0 and extend != 1:
            raise Exception("Wrong Input to extend parameter")
        
        
        blocked_mle_zero = self.blocked_reactions + self.mle_reactions + self.zero_flux_reactions
        list_removed_reactions = list(set(blocked_mle_zero))
        
        self.removed_reactions = list_removed_reactions
   
        
        self.remove_model_reactions()
        
                
        remained_reactions = list((Counter(self.initial_reactions)-Counter(self.removed_reactions)).elements())
        remained_reactions = list((Counter(remained_reactions)-Counter(self.essential_reactions)).elements())
   
        additional_removed_reactions_count = 0
        
        for reaction in remained_reactions:
            
            fba_solution_before = self.model.optimize().objective_value
                        
            initial_lower = self.model.reactions.get_by_id(reaction).lower_bound
            initial_upper = self.model.reactions.get_by_id(reaction).upper_bound
            
            self.model.reactions.get_by_id(reaction).lower_bound = 0
            self.model.reactions.get_by_id(reaction).upper_bound = 0
    
            fba_solution_after = self.model.optimize().objective_value
            
            if fba_solution_after != None and (extend == 1):
                if (abs(fba_solution_after - fba_solution_before) < tol):
                    self.removed_reactions.append(reaction)
                    additional_removed_reactions_count += 1
              
            self.model.reactions.get_by_id(reaction).upper_bound = initial_upper
            self.model.reactions.get_by_id(reaction).lower_bound = initial_lower
            
            
        fba_solution_initial = model.optimize().objective_value
        self.remove_model_reactions()        
        fba_solution_final = model.optimize().objective_value
        
        additional_removed_reactions_list = (self.removed_reactions[len(self.removed_reactions)-additional_removed_reactions_count:])
        
        if (fba_solution_final == None):
            for reaction in additional_removed_reactions_list:
                self.model.reactions.get_by_id(reaction).bounds = self.reaction_bounds_dict[reaction]
                self.removed_reactions.remove(reaction)
            print(len(self.removed_reactions), "of the", len(self.initial_reactions), "reactions were removed from the model")

        elif(abs(fba_solution_final - fba_solution_initial) > tol):
            for reaction in additional_removed_reactions_list:
                self.model.reactions.get_by_id(reaction).bounds = self.reaction_bounds_dict[reaction]
                self.removed_reactions.remove(reaction)
            print(len(self.removed_reactions), "of the", len(self.initial_reactions), "reactions were removed from the model") 

        else:
            print(len(self.removed_reactions), "of the", len(self.initial_reactions), "reactions were removed from the model")      
        
        
        return self.removed_reactions
     
        

model = load_json_model("ext_data/e_coli_core.json")
model = load_json_model("../../../iAF1260.json")

fba_solution = model.optimize()
print(fba_solution.objective_value)

obj = PreProcess(model)
rem = obj.removed(extend=1)
print(len(rem))

fba_solution = model.optimize()
print(fba_solution.objective_value)


