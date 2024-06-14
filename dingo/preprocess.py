
import cobra
from cobra.io import load_json_model

model = load_json_model("../ext_data/e_coli_core.json")


model.objective = 'BIOMASS_Ecoli_core_w_GAM'
fba_solution = model.optimize()
print(fba_solution.objective_value)


model.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM").lower_bound = fba_solution.objective_value
print(model.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM").lower_bound)
print(model.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM").upper_bound)


fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.95)
blocked = cobra.flux_analysis.find_blocked_reactions(model)
print(len(blocked))