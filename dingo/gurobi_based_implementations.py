# Dingo : a python library for metabolic networks sampling and analysis
# Dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB


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
                model.addMConstr(Aeq_sparse, x, '=', beq, name = "c")

                # Update the model to include the constraints added
                model.update()

                # Add constraints for the uneqalities of A
                model.addMConstr(A_sparse, x, '<', b, name = "d")

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





# Build a Python function to perform fva using gurobi LP solver
def fast_fva(lb, ub, S):

    if (lb.size != S.shape[1] or ub.size != S.shape[1]):
        raise Exception('The number of reactions must be equal to the number of given flux bounds.')

    m = S.shape[0]; n = S.shape[1]
    beq = np.zeros(m)

    A = np.zeros((2*n, n), dtype='float')
    A[0:n] = np.eye(n)
    A[n:] -=  np.eye(n,n, dtype='float')

    b = np.concatenate((ub, -lb), axis=0)
    b = np.asarray(b, dtype = 'float')
    b = np.ascontiguousarray(b, dtype = 'float')

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
                model.addMConstr(Aeq_sparse, x, '=', beq, name = "c")

                # Update the model to include the constraints added
                model.update()

                # Add constraints for the uneqalities of A
                model.addMConstr(A_sparse, x, '<', b, name = "d")

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


# This is a function to get the maximum ball inscribed in the polytope using gurobi solver
def fast_inner_ball(A, b):

   extra_column = []

   m = A.shape[0]
   n = A.shape[1]

   for i in range(A.shape[0]):
      entry = np.linalg.norm(A[i,])
      extra_column.append(entry)

   column = np.asarray(extra_column)
   A_expand = np.c_[A, column]

   with gp.Env(empty=True) as env:
      env.setParam('OutputFlag', 0)
      env.start()

      d = A_expand.shape[1]

      with gp.Model(env=env) as model:

         # Create variable
         x = model.addMVar(shape = d, vtype = GRB.CONTINUOUS , name = "x", lb = -GRB.INFINITY, ub = GRB.INFINITY)
         model.update()

         # Make A_full_dim sparse
         A_expand_sparse = sp.csr_matrix(A_expand.astype('float'))

         # Add constraints
         model.addMConstr(A_expand_sparse, x, '<', b, name = "c")
         model.update()

         # Set the ith row of the A matrix as the objective function
         a = np.ones((1,n+1))
         azero = np.zeros((1,n))
         a[:,:-1] = azero
         objective_function = a[0]

         # Set the objective function in the model
         model.setMObjective(None, objective_function, 0.0, None, None, x, GRB.MAXIMIZE)
         model.update()

         # Optimize model
         model.optimize()

         # Get the solution returned
         varss = model.getVars()

         # Get the center point and the radius of max ball from the solution of LP; its last element is the radius
         point = []
         for i in range(len(varss)):
            if i == len(varss) - 1:
               r = varss[i].x
            else:
               value = varss[i].x
               point.append(value)

         # And check whether the computed radius is negative
         if r < 0 :
            print ("The radius calculated has negative value. The polytope is infeasible or something went wrong with the solver")
         else:
            return point, r

