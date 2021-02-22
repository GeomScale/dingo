import numpy as np
from scipy.optimize import linprog

# For the fva, we need the following dependencies
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB


# This is a function to get the maximum ball inscribed in the polytope using scipy.optimize LP solver `linprog`
def slow_inner_ball(A, b):

    extra_column = []

    m = A.shape[0]
    n = A.shape[1]

    for i in range(A.shape[0]):
        entry = np.linalg.norm(A[i,])
        extra_column.append(entry)

    column = np.asarray(extra_column)
    A_expand = np.c_[A, column]
    
    a = np.ones((1,n+1))
    azero = np.zeros((1,n))
    a[:,:-1] = azero
    objective_function = a[0]

    res = linprog(objective_function, A_ub = A_expand, b_ub = b, bounds = (None, None))

    sol = res.x

    # Get the center point and the radius of max ball from the solution of LP; its last element is the radius
    point = []
    for i in range(len(sol)):
        if i == len(sol) - 1:
            r = sol[i].x
        else:
            value = sol[i].x
            point.append(value)
    
    # And check whether the computed radius is negative
    if r < 0 :
        print ("The radius calculated has negative value. The polytope is infeasible or something went wrong with the solver")
    else:
        return point, r



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
   print("A_expand :") ; print(A_expand)
   print("A_expand shape :") ; print(A_expand.shape)

   with gp.Env(empty=True) as env:
      env.setParam('OutputFlag', 0)
      env.start()

      d = A_expand.shape[1]

      with gp.Model(env=env) as model:

         # Create variable
         x = model.addMVar(shape = d, vtype = GRB.CONTINUOUS , name = "x", lb = -GRB.INFINITY, ub = GRB.INFINITY)
         model.update()

         # Make A_full_dim sparse
         A_expand_sparse = sp.csr_matrix(A_expand.astype(np.float))

         # Add constraints
         model.addMConstrs(A_expand_sparse, x, '<', b, name = "c")
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
  
  
