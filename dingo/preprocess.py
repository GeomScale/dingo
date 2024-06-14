
import cobra
from cobra.io import load_json_model


class PreProcess:
    
    def __init__(self, model):
        self.model = model
    
    def blocked_reactions(model):
        
        tol = 1e-6

        model.objective = 'BIOMASS_Ecoli_core_w_GAM'
        fba_solution = model.optimize()

        model.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM").lower_bound = fba_solution.objective_value

        fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.95)
        blocked_fva = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        blocked_reactions = blocked_fva.index.tolist()
        
        return blocked_reactions


model = load_json_model("../ext_data/e_coli_core.json")
blocked = PreProcess.blocked_reactions(model)
print(blocked)

