
import cobra
from cobra.io import load_json_model


class PreProcess:
    
    def __init__(self, model):
        self.model = model
        
    def objective_function(model):
        
        objective = str(model.summary()._objective)
        objective = objective.split(" ")[1]
        return objective

    def metabolically_less_efficient(model):
        
        objective = PreProcess.objective_function(model)
    
        tol = 1e-6

        model.objective = objective
        fba_solution = model.optimize()

        model.reactions.get_by_id(objective).lower_bound = fba_solution.objective_value

        fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.95)
        blocked_fva = fva.loc[ (abs(fva['minimum']) < tol ) & (abs(fva['maximum']) < tol)]
        mle = blocked_fva.index.tolist()
        
        return mle


model = load_json_model("../ext_data/e_coli_core.json")

blocked = PreProcess.metabolically_less_efficient(model)
print(blocked)

