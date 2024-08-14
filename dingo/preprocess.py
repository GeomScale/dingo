
import cobra
import cobra.manipulation
from collections import Counter
from dingo import MetabolicNetwork, PolytopeSampler
from dingo.utils import correlated_reactions
import numpy as np


class PreProcess:
    
    def __init__(self, model, tol = 1e-6, open_exchanges = False, verbose = False):

        """
        model -- parameter gets a cobra model as input
        
        tol -- parameter gets a cutoff value used to classify 
               zero-flux and mle reactions and compare FBA solutions
               before and after reactions removal
        
        open_exchanges -- parameter is used in the function that identifies blocked reactions
                          It controls whether or not to open all exchange reactions 
                          to very high flux ranges.
                          
        verbose -- A boolean type variable that if True 
                   additional information for preprocess is printed.
        """
        
        self._model = model
        self._tol = tol
                
        self._open_exchanges = open_exchanges
        self._verbose = verbose
        
        if self._tol > 1e-6 and verbose == True:
            print("Tolerance value set to",self._tol,"while default value is 1e-6. A looser check will be performed")

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
        self._objective = objective.split(" ")[1]
        
        return self._objective


    def _initial(self):
        """
        A function used to find reaction ids of a model
        """
        
        self._initial_reactions =  [ reaction.id for reaction in \
        self._model.reactions ]
                    
        return self._initial_reactions
    
    
    def _reaction_bounds_dictionary(self):
        """
        A function used to create a dictionary that maps
        reactions with their corresponding bounds. It is used to
        later restore bounds to their wild-type values
        """
        
        self._reaction_bounds_dict = { }
                
        for reaction_id in self._initial_reactions:
            bounds = self._model.reactions.get_by_id(reaction_id).bounds
            self._reaction_bounds_dict[reaction_id] = bounds
            
        return self._reaction_bounds_dict


    def _essentials(self):
        """
        A function used to find all the essential reactions
        and append them into a list. Essential reactions are
        the ones that are required for growth. If removed the 
        objective function gets zeroed.
        """
        
        self._essential_reactions =  [ reaction.id for reaction in \
        cobra.flux_analysis.find_essential_reactions(self._model) ]
                    
        return self._essential_reactions
    

    def _zero_flux(self):
        """
        A function used to find zero-flux reactions.
        “Zero-flux” reactions cannot carry a flux while maintaining 
        at least 90% of the maximum growth rate.
        These reactions have both a min and a max flux equal to 0,
        when running a FVA analysis with the fraction of optimum set to 90%
        """
        
        tol = self._tol
        
        fva = cobra.flux_analysis.flux_variability_analysis(self._model, fraction_of_optimum=0.9)
        zero_flux = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        self._zero_flux_reactions = zero_flux.index.tolist()
        
        return self._zero_flux_reactions
    
    
    def _blocked(self):
        """
        A function used to find blocked reactions.
        "Blocked" reactions that cannot carry a flux in any condition.
        These reactions can not have any flux other than 0
        """
        
        self._blocked_reactions = cobra.flux_analysis.find_blocked_reactions(self._model, open_exchanges=self._open_exchanges)
        return self._blocked_reactions


    def _metabolically_less_efficient(self):
        """
        A function used to find metabolically less efficient reactions.
        "Metabolically less efficient" require a reduction in growth rate if used
        These reactions are found when running an FBA and setting the  
        optimal growth rate as the lower bound of the objective function (in 
        this case biomass production. After running an FVA with the fraction of optimum
        set to 0.95, the reactions that have no flux are the metabolically less efficient.
        """
            
        tol = self._tol

        fba_solution = self._model.optimize()

        wt_lower_bound = self._model.reactions.get_by_id(self._objective).lower_bound
        self._model.reactions.get_by_id(self._objective).lower_bound = fba_solution.objective_value

        fva = cobra.flux_analysis.flux_variability_analysis(self._model, fraction_of_optimum=0.95)
        mle = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        self._mle_reactions = mle.index.tolist()
        
        self._model.reactions.get_by_id(self._objective).lower_bound = wt_lower_bound
        
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
   
            
    def reduce(self, extend=False):
        """
        A function that calls the "remove_model_reactions" function
        and removes blocked, zero-flux and metabolically less efficient 
        reactions from the model.
        
        Then it finds the remaining reactions in the model after 
        exclusion of the essential reactions.
        
        When the "extend" parameter is set to True, the function  performes
        an additional check to remove further reactions. These reactions 
        are the ones that if knocked-down, they do not affect the value 
        of the objective function. Reactions are removed in an ordered way. 
        The ones with the least overall correlation in a correlation matrix 
        are removed first. These reactions are removed one by one from the model. 
        If this removal produces an infeasible solution (or a solution of 0) 
        to the objective function, these reactions are restored to their initial bounds.
        
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
        
        tol = self._tol
        
        if extend != False and extend != True:
            raise Exception("Wrong Input to extend parameter")
        
        elif extend == False:
            
            if self._verbose == True:
                print(len(self._removed_reactions), "of the", len(self._initial_reactions), \
                "reactions were removed from the model with extend set to", extend) 
            
            # call this functon to convert cobra to dingo model
            self._dingo_model = MetabolicNetwork.from_cobra_model(self._model)
            return self._removed_reactions, self._dingo_model 

        elif extend == True:
            
            reduced_dingo_model = MetabolicNetwork.from_cobra_model(self._model)
            reactions = reduced_dingo_model.reactions      
            sampler = PolytopeSampler(reduced_dingo_model)
            steady_states = sampler.generate_steady_states()

            # calculate correlation matrix with additional filtering from copula indicator       
            corr_matrix = correlated_reactions(
                                steady_states,  
                                pearson_cutoff = 0,
                                indicator_cutoff = 0, 
                                cells = 10,
                                cop_coeff = 0.3,
                                lower_triangle = False)
            
            # convert pearson values to absolute values
            abs_array = abs(corr_matrix)
            # sum absolute pearson values per row 
            sum_array = np.sum((abs_array), axis=1)
            # get indices of ordered sum values
            order_sum_indices = np.argsort(sum_array)
            
            fba_solution_before = self._model.optimize().objective_value
            
            # count additional reactions with a possibility of removal
            additional_removed_reactions_count = 0 
            
            # find additional reactions with a possibility of removal
            for index in order_sum_indices:
                if reactions[index] in remained_reactions:
                    reaction = reactions[index]
                                
                    # perform a knock-out and check the output
                    self._model.reactions.get_by_id(reaction).lower_bound = 0
                    self._model.reactions.get_by_id(reaction).upper_bound = 0
    
                    try:
                        fba_solution_after = self._model.optimize().objective_value
                        if (abs(fba_solution_after - fba_solution_before) > tol):
                            # restore bounds
                            self._model.reactions.get_by_id(reaction).bounds = self._reaction_bounds_dict[reaction]
                    
                    # if system has no solution
                    except:
                        self._model.reactions.get_by_id(reaction).bounds = self._reaction_bounds_dict[reaction]
                    
                    finally:
                        if (fba_solution_after != None) and (abs(fba_solution_after - fba_solution_before) < tol):
                            self._removed_reactions.append(reaction)
                            additional_removed_reactions_count += 1
                        else:
                            # restore bounds
                            self._model.reactions.get_by_id(reaction).bounds = self._reaction_bounds_dict[reaction]

            
            if self._verbose == True:  
                print(len(self._removed_reactions), "of the", len(self._initial_reactions), \
                "reactions were removed from the model with extend set to", extend)
                print(additional_removed_reactions_count, "additional reaction(s) removed")
                        
            # call this functon to convert cobra to dingo model
            self._dingo_model = MetabolicNetwork.from_cobra_model(self._model)
            return self._removed_reactions, self._dingo_model 
