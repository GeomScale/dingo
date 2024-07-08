
import cobra
import cobra.manipulation
from collections import Counter
from dingo.loading_models import parse_cobra_model


class PreProcess:
    
    def __init__(self, model):
        self.model = model
        self.objective = self.objective_function()
        self.initial_reactions = self.initial()
        self.reaction_bounds_dict = self.reaction_bounds_dictionary()
        self.essential_reactions = self.essentials()
        self.zero_flux_reactions = self.zero_flux()
        self.blocked_reactions = self.blocked()
        self.mle_reactions = self.metabolically_less_efficient()
        self.dingo_model = self.cobra_dingo_conversion()
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
        """
        A function used to find reaction ids of a model
        """
        
        self.initial_reactions = []
        
        for reaction in self.model.reactions:
            reaction_id = reaction.id
            self.initial_reactions.append(reaction_id)
            
        return self.initial_reactions
    
    
    def reaction_bounds_dictionary(self):
        """
        A function used to create a dictionary that maps
        reactions with reactions bounds. It is used to
        later restore some bounds to wild-type values
        """
        
        self.reaction_bounds_dict = {
                  'reaction': (0, 100)
                  }
                
        for reaction_id in self.initial_reactions:
            bounds = self.model.reactions.get_by_id(reaction_id).bounds
            self.reaction_bounds_dict[reaction_id] = bounds
            
        return self.reaction_bounds_dict


    def essentials(self):
        """
        A function used to find all essential reactions
        and appends them in a list
        """
        
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
        """
        A function used to find metabolically less efficient reactions.
        These reactions are found when running an FBA and setting the  
        optimal growth rate as the lower bound of the objective function (in 
        this case biomass production. After running an FVA with a fraction of optimum
        set to 0.95, the reactions that have no flux are the metabolically less efficient.
        """
            
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
        """
        A function used to set lower and upper bounds of certain reactions to 0
        (it turns off reactions)
        """
                                
        for reaction in self.removed_reactions:
            self.model.reactions.get_by_id(reaction).lower_bound = 0
            self.model.reactions.get_by_id(reaction).upper_bound = 0
            
        return self.model
    
    
    def cobra_dingo_conversion(self):
        """
        A function used to convert the reduced cobra model to a dingo model
        """
        self.dingo_model = parse_cobra_model(self.model)
        return self.dingo_model
   
            
    def reduce(self, extend):
        """
        A function that calls "remove_model_reactions" function
        and removes blocked, zero-flux and metabolically less efficient reactions.
        Then it finds the remaining reactions in the model after 
        exclusion of the essential reactions.
        
        The "extend" parameter when set to 1 performes an additional check to remove
        further reactions. These reactions are the ones that if knocked-down, they
        do not affect the value of the objective function. These reactions 
        are removed simultaneously. If the simultaneous removal produces an infesible
        solution (or 0) to the objective function, they are restored with their initial bounds.
        
        The dingo model is then created from the cobra model 
        using the "cobra_dingo_conversion" function.
        
        The outputs are
        (a) A list of the removed reactions ids
        (b) The reduced dingo model
        """        
        
        # create a list from the combined blocked, zero-flux, mle reactions
        blocked_mle_zero = self.blocked_reactions + self.mle_reactions + self.zero_flux_reactions
        list_removed_reactions = list(set(blocked_mle_zero))        
        self.removed_reactions = list_removed_reactions
   
        # remove these reactions from the model
        self.remove_model_reactions()
                     
        remained_reactions = list((Counter(self.initial_reactions)-Counter(self.removed_reactions)).elements())
        remained_reactions = list((Counter(remained_reactions)-Counter(self.essential_reactions)).elements())
        
        tol = 1e-6
        
        if extend != 0 and extend != 1:
            raise Exception("Wrong Input to extend parameter")
  
        # find additional reactions with a possibility of removal
        additional_removed_reactions_count = 0       
        for reaction in remained_reactions:
            
            fba_solution_before = self.model.optimize().objective_value
                        
            initial_lower = self.model.reactions.get_by_id(reaction).lower_bound
            initial_upper = self.model.reactions.get_by_id(reaction).upper_bound
            
            # perform a knock-out and check the output
            self.model.reactions.get_by_id(reaction).lower_bound = 0
            self.model.reactions.get_by_id(reaction).upper_bound = 0
    
            fba_solution_after = self.model.optimize().objective_value
            
            if fba_solution_after != None and (extend == 1):
                if (abs(fba_solution_after - fba_solution_before) < tol):
                    self.removed_reactions.append(reaction)
                    additional_removed_reactions_count += 1
              
            self.model.reactions.get_by_id(reaction).upper_bound = initial_upper
            self.model.reactions.get_by_id(reaction).lower_bound = initial_lower
            
            
        # compare FBA solution before and after the removal of additional reactions
        fba_solution_initial = self.model.optimize().objective_value
        self.remove_model_reactions()        
        fba_solution_final = self.model.optimize().objective_value

        
        additional_removed_reactions_list = (self.removed_reactions[len(self.removed_reactions)-additional_removed_reactions_count:])
        
        # if FBA solution after removal is infesible or altered
        # restore the initial reactions bounds
        if (fba_solution_final == None):
            for reaction in additional_removed_reactions_list:
                self.model.reactions.get_by_id(reaction).bounds = self.reaction_bounds_dict[reaction]
                self.removed_reactions.remove(reaction)
            print(len(self.removed_reactions), "of the", len(self.initial_reactions), "reactions were removed from the model with extend set to", extend)

        elif(abs(fba_solution_final - fba_solution_initial) > tol):
            for reaction in additional_removed_reactions_list:
                self.model.reactions.get_by_id(reaction).bounds = self.reaction_bounds_dict[reaction]
                self.removed_reactions.remove(reaction)
            print(len(self.removed_reactions), "of the", len(self.initial_reactions), "reactions were removed from the model with extend set to", extend) 

        else:
            print(len(self.removed_reactions), "of the", len(self.initial_reactions), "reactions were removed from the model with extend set to", extend)
            
        
        return self.removed_reactions, self.cobra_dingo_conversion() 
