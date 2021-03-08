import numpy as np
from scipy.optimize import linprog

# For the fva, we need the following dependencies
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB


import sys
# Build a Python function to perform fva using scipy.optimize LP solver `linprog`
def slow_fva(lb, ub, S):

    if (lb.size != S.shape[1] or ub.size != S.shape[1]):
        raise Exception('The number of reactions must be equal to the number of given flux bounds.')
    if (c.size != S.shape[1]):
        raise Exception('The length of the lineart objective function must be equal to the number of reactions.')

    m = S.shape[0]; n = S.shape[1]
    optimum_value = 0
    optimum_sol = np.zeros(n)

    A = np.zeros((2*n, n), dtype=np.float)
    A[0:n] = np.eye(n)
    A[n:] -= np.eye(n, n, dtype=np.float)

    b = np.concatenate((lb, -ub), axis=1)
    beq = np.zeros(m)

    Aeq_new = S
    beq_new = beq
   
    min_fluxes = []
    max_fluxes = []

    try:

        for i in range(int(A.shape[0]/2)):
            
            # Set the ith row of the A matrix as the objective function
            objective_function = A[i,]
            res = linprog(objective_function, A_ub = A, b_ub = b, A_eq = S, b_eq = beq, bounds = (None, None))

            # If optimized
            if res.success:
                min_objective = res.fun
                min_fluxes.append(min_objective)
                
            
            # Set the minus of the ith row of the A matrix as the objective function
            objective_function = np.asarray([-x for x in objective_function])
            res = linprog(objective_function, A_ub = A, b_ub = b, A_eq = S, b_eq = beq, bounds = (None, None))

            # If optimized
            if res.success:
                max_objective = -res.fun
                max_fluxes.append(max_objective)
                
            
            # Calculate the width
            width = abs(max_objective - min_objective)                  

            # Check whether we keep or not the equality
            if width < 1e-07:
                Aeq_new = np.vstack((Aeq_new, A[i,]))
                beq_new = np.append(beq_new, min(max_objective, min_objective))
        
        # The np.vstack() creates issues on changing contiguous c orded of np arrays; here we fix this
        Aeq_new = np.ascontiguousarray(Aeq_new, dtype=np.dtype)

        # Furthremore, we need to have float64 in all numpy arrays
        Aeq_new = Aeq_new.astype('float64')

        # Make lists of fluxes numpy arrays
        min_fluxes = np.asarray(min_fluxes) ; max_fluxes = np.asarray(max_fluxes)

        return A, b, Aeq_new, beq_new, min_fluxes, max_fluxes

    except AttributeError :
        print ("scipy.optimize.linprog failed.")



# Build a Python function to perform fva using gurobi LP solver
def fast_fva(lb, ub, S):

    if (lb.size != S.shape[1] or ub.size != S.shape[1]):
        raise Exception('The number of reactions must be equal to the number of given flux bounds.')
    if (c.size != S.shape[1]):
        raise Exception('The length of the lineart objective function must be equal to the number of reactions.')

    m = S.shape[0]; n = S.shape[1]
    Aeq_new = S
    beq_new = beq
   
    min_fluxes = []
    max_fluxes = []

    try:

        # To avoid printint the output of the optimize() function of Gurobi, we need to set an environment like this
        with gp.Env(empty=True) as env:
            env.setParam('OutputFlag', 0)
            env.start()

            with gp.Model(env=env) as model:

                # Create variables
                x = model.addMVar(shape = n, vtype = GRB.CONTINUOUS , name = "x", lb = -GRB.INFINITY, ub = GRB.INFINITY)

                # Make sparse Aeq
                Aeq_sparse = sp.csr_matrix(S)

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

                # Loop through the lines of the A matrix, set objective function for each and run the model
                for i in range(int(A.shape[0]/2)):
            
                    # Set the ith row of the A matrix as the objective function
                    objective_function = A[i,]

                    # Set the objective function in the model
                    model.setMObjective(None, objective_function, 0.0, None, None, x, GRB.MINIMIZE)
                    model.update()

                    # Optimize model
                    model.optimize()

                    # If optimized
                    status = model.status
                    if status == GRB.OPTIMAL:

                        # Get the min objective value
                        min_objective = model.getObjective().getValue()
                        min_fluxes.append(min_objective)

                    # Likewise, for the maximum
                    objective_function = np.asarray([-x for x in objective_function])
                    model.setMObjective(None, objective_function, 0.0, None, None, x, GRB.MINIMIZE)
                    model.update()
                    model.optimize()

                    # Again if optimized
                    status = model.status
                    if status == GRB.OPTIMAL:

                        # Get the max objective value
                        max_objective = -model.getObjective().getValue()
                        max_fluxes.append(max_objective)

                    # Calculate the width
                    width = abs(max_objective - min_objective)                  

                    # Check whether we keep or not the equality
                    if width < 1e-07:
                        Aeq_new = np.vstack((Aeq_new, A[i,]))
                        beq_new = np.append(beq_new, min(max_objective, min_objective))                                  
               
                # The np.vstack() creates issues on changing contiguous c orded of np arrays; here we fix this
                Aeq_new = np.ascontiguousarray(Aeq_new, dtype=np.dtype)

                # Furthremore, we need to have float64 in all numpy arrays
                Aeq_new = Aeq_new.astype('float64')

                # Make lists of fluxes numpy arrays
                min_fluxes = np.asarray(min_fluxes) ; max_fluxes = np.asarray(max_fluxes)

                return A, b, Aeq_new, beq_new, min_fluxes, max_fluxes


    # Print error messages
    except gp.GurobiError as e :
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError :
        print ("Gurobi solver failed.")
  
  
  
