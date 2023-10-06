import gurobipy as gp
import numpy as np
import scipy.sparse as sp

def update_model_constraints_and_bounds(model, Aeq_sparse=None, beq=None, A_sparse=None, b=None,
                                       new_lower_bounds=None, new_upper_bounds=None):
    if Aeq_sparse is not None and beq is not None:
        # Remove the old equality constraints
        model.remove("c")

        # Add new equality constraints
        model.addMConstr(sp.csr_matrix(Aeq_sparse), model.getVars(), "=", np.array(beq), name="c")

    if A_sparse is not None and b is not None:
        # Remove the old inequality constraints
        model.remove("d")

        # Add new inequality constraints
        model.addMConstr(sp.csr_matrix(A_sparse), model.getVars(), "<", np.array(b), name="d")

    if new_lower_bounds is not None or new_upper_bounds is not None:
        # Retrieve the variable objects
        x_vars = model.getVars()

        # Update variable bounds if specified
        if new_lower_bounds is not None:
            for i, x_var in enumerate(x_vars):
                x_var.lb = new_lower_bounds[i]

        if new_upper_bounds is not None:
            for i, x_var in enumerate(x_vars):
                x_var.ub = new_upper_bounds[i]

    # Update the model to apply the changes
    model.update()

# Usage example:
# Assuming you have a Gurobi model 'model' and new constraint matrices/vectors and bounds
# Call the function to update the model
# update_model_constraints_and_bounds(model, Aeq_sparse_new, beq_new, A_sparse_new, b_new, new_lower_bounds, new_upper_bounds)

