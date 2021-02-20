import numpy as np
from scipy.optimize import linprog

# For the fva, we need the following dependencies
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB



# Build a Python function to perform fba using scipy.optimize LP solver `linprog`
def slow_fba(A, b, Aeq, beq, c):

    d = A.shape[1] ; m = Aeq.shape[0] ; n = Aeq.shape[1]
    optimum_value = 0
    optimum_sol = np.zeros(d)

    try:

        #objective_function = np.asarray(c)
        objective_function = np.asarray([-x for x in c])

        res = linprog(objective_function, A_ub = A, b_ub = b, A_eq = Aeq, b_eq = beq, bounds = (None, None))

        # If optimized
        if res.success:
            optimum_value = -res.fun
            optimum_sol = np.asarray(res.x)

        return optimum_sol, optimum_value

    except AttributeError :
        print ("Encountered an attribute error ")




# Build a Python function to perform fba using gurobi LP solver
def fast_fba(A, b, Aeq, beq, c):

    d = A.shape[1] ; m = Aeq.shape[0] ; n = Aeq.shape[1]
    optimum_value = 0
    optimum_sol = np.zeros(d)

    try:

      # To avoid printint the output of the optimize() function of Gurobi, we need to set an environment like this
        with gp.Env(empty=True) as env:
            env.setParam('OutputFlag', 0)
            env.start()

            with gp.Model(env=env) as model:

                # Create variables
                x = model.addMVar(shape = d, vtype = GRB.CONTINUOUS , name = "x", lb = -GRB.INFINITY, ub = GRB.INFINITY)

                # Make sparse Aeq
                Aeq_sparse = sp.csr_matrix(Aeq)

                # Make A sparse
                A_sparse = sp.csr_matrix(A)            

                # Set the b and beq vectors as numpy vectors
                b = np.array(b)
                beq = np.array(beq)

                # Add constraints
                model.addMConstrs(Aeq_sparse, x, '=', beq, name = "c")

                # Update the model to include the constraints added
                model.update()

                # Add constraints for the uneqalities of A
                model.addMConstrs(A_sparse, x, '<', b, name = "d")

                # Update the model with the extra constraints and then print it
                model.update()

                objective_function = np.asarray(c)
                model.setMObjective(None, objective_function, 0.0, None, None, x, GRB.MAXIMIZE)

                # If optimized
                status = model.status
                if status == GRB.OPTIMAL:
                    optimum_value = -model.getObjective().getValue()
                    optimum_sol = np.asarray(model.getObjective())

                return optimum_sol, optimum_value


        # Print error messages
    except gp . GurobiError as e :
        print ("Error code " + str( e . errno ) + ": " + str( e ))
    except AttributeError :
        print ("Encountered an attribute error ")
  
  
