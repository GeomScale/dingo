
import cobra
import cobra.manipulation
from collections import Counter
from dingo import MetabolicNetwork


class PreProcess:
    
    def __init__(self, model):

        self._model = model
        self._objective = self._objective_function()
        self._initial_reactions = self._initial()
        self._reaction_bounds_dict = self._reaction_bounds_dictionary()
        self._essential_reactions = self._essentials()
        self._zero_flux_reactions = self._zero_flux()
        self._blocked_reactions = self._blocked()
        self._mle_reactions = self._metabolically_less_efficient()
        self._removed_reactions = []


    def _objective_function(self):
        """
        A function used to find the objective function of a model
        """
        
        objective = str(self._model.summary()._objective)
        objective = objective.split(" ")[1]
        
        self._objective = objective
        return self._objective


    def _initial(self):
        """
        A function used to find reaction ids of a model
        """
        
        self._initial_reactions = []
        
        for reaction in self._model.reactions:
            reaction_id = reaction.id
            self._initial_reactions.append(reaction_id)
            
        return self._initial_reactions
    
    
    def _reaction_bounds_dictionary(self):
        """
        A function used to create a dictionary that maps
        reactions with their corresponding bounds. It is used to
        later restore bounds to their wild-type values
        """
        
        self._reaction_bounds_dict = {
                  'reaction': (0, 1000)
                  }
                
        for reaction_id in self._initial_reactions:
            bounds = self._model.reactions.get_by_id(reaction_id).bounds
            self._reaction_bounds_dict[reaction_id] = bounds
            
        return self._reaction_bounds_dict


    def _essentials(self):
        """
        A function used to find all the essential reactions
        and append them into a list
        """
        
        self._essential_reactions = []
        essential_reactions = cobra.flux_analysis.find_essential_reactions(self._model)
        for reaction in essential_reactions:
            reaction_id = reaction.id
            self._essential_reactions.append(reaction_id)
            
        return self._essential_reactions
    

    def _zero_flux(self):
        """
        A function used to find zero-flux reactions.
        These reactions are the ones that have a flux equaled to 0
        when running a FVA analysis with the fraction of optimum set to 90%
        """
        
        tol = 1e-6
        
        fva = cobra.flux_analysis.flux_variability_analysis(self._model, fraction_of_optimum=0.9)
        zero_flux = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        zero_flux_reactions = zero_flux.index.tolist()
        
        self._zero_flux_reactions = zero_flux_reactions
        return self._zero_flux_reactions
    
    
    def _blocked(self):
        """
        A function used to find blocked reactions.
        These reactions can not have any flux other than 0
        """
        
        blocked_reactions = cobra.flux_analysis.find_blocked_reactions(self._model)
        self._blocked_reactions =  blocked_reactions
        return self._blocked_reactions


    def _metabolically_less_efficient(self):
        """
        A function used to find metabolically less efficient reactions.
        These reactions are found when running an FBA and setting the  
        optimal growth rate as the lower bound of the objective function (in 
        this case biomass production. After running an FVA with the fraction of optimum
        set to 0.95, the reactions that have no flux are the metabolically less efficient.
        """
            
        tol = 1e-6

        self._model.objective = self._objective
        fba_solution = self._model.optimize()

        wt_lower_bound = self._model.reactions.get_by_id(self._objective).lower_bound
        self._model.reactions.get_by_id(self._objective).lower_bound = fba_solution.objective_value

        fva = cobra.flux_analysis.flux_variability_analysis(self._model, fraction_of_optimum=0.95)
        mle = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        mle = mle.index.tolist()
        
        self._model.reactions.get_by_id(self._objective).lower_bound = wt_lower_bound
        
        self._mle_reactions = mle
        return self._mle_reactions
    
    
    def _remove_model_reactions(self):
        """
        A function used to set the lower and upper bounds of certain reactions to 0
        (it turns off reactions)
        """
                                
        for reaction in self._removed_reactions:
            self._model.reactions.get_by_id(reaction).lower_bound = 0
            self._model.reactions.get_by_id(reaction).upper_bound = 0
            
        return self._model
   
            
    def reduce(self, extend):
        """
        A function that calls the "remove_model_reactions" function
        and removes blocked, zero-flux and metabolically less efficient 
        reactions from the model.
        
        Then it finds the remaining reactions in the model after 
        exclusion of the essential reactions.
        
        When the "extend" parameter is set to 1, the function  performes
        an additional check to remove further reactions. These reactions 
        are the ones that if knocked-down, they do not affect the value 
        of the objective function. These reactions are removed simultaneously
        from the model. If this simultaneous removal produces an infesible
        solution (or a solution of 0) to the objective function, 
        these reactions are restored with their initial bounds.
        
        A dingo-type tuple is then created from the cobra model 
        using the "cobra_dingo_tuple" function.
        
        The outputs are
        (a) A list of the removed reactions ids
        (b) A reduced dingo model
        """        
        
        # create a list from the combined blocked, zero-flux, mle reactions
        blocked_mle_zero = self._blocked_reactions + self._mle_reactions + self._zero_flux_reactions
        list_removed_reactions = list(set(blocked_mle_zero))        
        self._removed_reactions = list_removed_reactions
   
        # remove these reactions from the model
        self._remove_model_reactions()
                     
        remained_reactions = list((Counter(self._initial_reactions)-Counter(self._removed_reactions)).elements())
        remained_reactions = list((Counter(remained_reactions)-Counter(self._essential_reactions)).elements())
        
        tol = 1e-6
        
        if extend != 0 and extend != 1:
            raise Exception("Wrong Input to extend parameter")
  
        # find additional reactions with a possibility of removal
        additional_removed_reactions_count = 0       
        for reaction in remained_reactions:
            
            fba_solution_before = self._model.optimize().objective_value
                        
            initial_lower = self._model.reactions.get_by_id(reaction).lower_bound
            initial_upper = self._model.reactions.get_by_id(reaction).upper_bound
            
            # perform a knock-out and check the output
            self._model.reactions.get_by_id(reaction).lower_bound = 0
            self._model.reactions.get_by_id(reaction).upper_bound = 0
    
            fba_solution_after = self._model.optimize().objective_value
            
            if fba_solution_after != None and (extend == 1):
                if (abs(fba_solution_after - fba_solution_before) < tol):
                    self._removed_reactions.append(reaction)
                    additional_removed_reactions_count += 1
              
            self._model.reactions.get_by_id(reaction).upper_bound = initial_upper
            self._model.reactions.get_by_id(reaction).lower_bound = initial_lower
            
            
        # compare FBA solution before and after the removal of additional reactions
        fba_solution_initial = self._model.optimize().objective_value
        self._remove_model_reactions()        
        fba_solution_final = self._model.optimize().objective_value

        
        # get a list of the additional reactions to be removed
        additional_removed_reactions_list = (self._removed_reactions[len(self._removed_reactions)-additional_removed_reactions_count:])
        
        # if FBA solution after removal is infesible or altered
        # restore the initial reactions bounds
        if (fba_solution_final == None):
            for reaction in additional_removed_reactions_list:
                self._model.reactions.get_by_id(reaction).bounds = self._reaction_bounds_dict[reaction]
                self._removed_reactions.remove(reaction)
            print(len(self._removed_reactions), "of the", len(self._initial_reactions), "reactions were removed from the model with extend set to", extend)

        elif(abs(fba_solution_final - fba_solution_initial) > tol):
            for reaction in additional_removed_reactions_list:
                self._model.reactions.get_by_id(reaction).bounds = self._reaction_bounds_dict[reaction]
                self.removed_reactions.remove(reaction)
            print(len(self._removed_reactions), "of the", len(self._initial_reactions), "reactions were removed from the model with extend set to", extend) 

        else:
            print(len(self._removed_reactions), "of the", len(self._initial_reactions), "reactions were removed from the model with extend set to", extend)

    
        # call this functon to convert cobra to dingo model
        self._dingo_model = MetabolicNetwork.from_cobra_model(self._model)

        return self._removed_reactions, self._dingo_model 
