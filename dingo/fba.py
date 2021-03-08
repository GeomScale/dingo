import numpy as np
from scipy.optimize import linprog

# For the fva, we need the following dependencies
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB


import sys
# Build a Python function to perform fba using scipy.optimize LP solver `linprog`
def slow_fba(lb, ub, S, c):

    if (lb.size != S.shape[1] or ub.size != S.shape[1]):
        raise Exception('The number of reactions must be equal to the number of given flux bounds.')
    if (c.size != S.shape[1]):
        raise Exception('The length of the lineart objective function must be equal to the number of reactions.')

    m = S.shape[0] ; n = S.shape[1]
    optimum_value = 0
    optimum_sol = np.zeros(n)

    A = np.zeros((2*n, n), dtype='float')
    A[0:n] = np.eye(n)
    A[n:] -= np.eye(n, n, dtype='float')

    b = np.concatenate((ub, -lb), axis=0)
    b = np.asarray(b, dtype = 'float')
    b = np.ascontiguousarray(b, dtype = 'float')

    beq = np.zeros(m)

    try:

        objective_function = np.asarray(c)
        objective_function = np.asarray([-x for x in c])

        res = linprog(objective_function, A_ub = A, b_ub = b, A_eq = S, b_eq = beq, bounds = (None, None))

        # If optimized
        if res.success:
            optimum_value = -res.fun
            optimum_sol = np.asarray(res.x)

        return optimum_sol, optimum_value

    except AttributeError :
        print ("scipy.optimize.linprog failed.")




# Build a Python function to perform fba using gurobi LP solver
def fast_fba(lb, ub, S, c):

    if (lb.size != S.shape[1] or ub.size != S.shape[1]):
        raise Exception('The number of reactions must be equal to the number of given flux bounds.')
    if (c.size != S.shape[1]):
        raise Exception('The length of the lineart objective function must be equal to the number of reactions.')

    m = S.shape[0]; n = S.shape[1]
    optimum_value = 0
    optimum_sol = []

    beq = np.zeros(m)

    A = np.zeros((2*n, n), dtype = 'float')
    A[0:n] = np.eye(n)
    A[n:] -=  np.eye(n,n, dtype = 'float')

    b = np.concatenate((ub, -lb), axis=0)
    b = np.asarray(b, dtype = 'float')
    b = np.ascontiguousarray(b, dtype = 'float')

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

                objective_function = np.asarray([-x for x in c])
                
                # Set the objective function in the model
                model.setMObjective(None, objective_function, 0.0, None, None, x, GRB.MINIMIZE)
                model.update()

                # Optimize model
                model.optimize()

                # If optimized
                status = model.status
                if status == GRB.OPTIMAL:
                    optimum_value = -model.getObjective().getValue()
                    v = model.getVars()

                for i in range(n):
                    optimum_sol.append(v[i].x)
                
                optimum_sol = np.asarray(optimum_sol)

                return optimum_sol, optimum_value


    # Print error messages
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    except AttributeError :
        print ("Gurobi failed.")
  
  
