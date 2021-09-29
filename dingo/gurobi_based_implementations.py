# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import numpy as np
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB
import math


def fast_fba(lb, ub, S, c):
    """A Python function to perform fba using gurobi LP solver
    Returns an optimal solution and its value for the following linear program:
    max c*v, subject to,
    Sv = 0, lb <= v <= ub

    Keyword arguments:
    lb -- lower bounds for the fluxes, i.e., a n-dimensional vector
    ub -- upper bounds for the fluxes, i.e., a n-dimensional vector
    S -- the mxn stoichiometric matrix, s.t. Sv = 0
    c -- the linear objective function, i.e., a n-dimensional vector
    """

    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )
    if c.size != S.shape[1]:
        raise Exception(
            "The length of the lineart objective function must be equal to the number of reactions."
        )
    tol = 1e-06

    m = S.shape[0]
    n = S.shape[1]
    optimum_value = 0
    optimum_sol = []

    beq = np.zeros(m)

    A = np.zeros((2 * n, n), dtype="float")
    A[0:n] = np.eye(n)
    A[n:] -= np.eye(n, n, dtype="float")

    b = np.concatenate((ub, -lb), axis=0)
    b = np.asarray(b, dtype="float")
    b = np.ascontiguousarray(b, dtype="float")

    try:

        # To avoid printint the output of the optimize() function of Gurobi, we need to set an environment like this
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()

            with gp.Model(env=env) as model:

                # Create variables
                x = model.addMVar(
                    shape=n,
                    vtype=GRB.CONTINUOUS,
                    name="x",
                    lb=-GRB.INFINITY,
                    ub=GRB.INFINITY,
                )

                # Make sparse Aeq
                Aeq_sparse = sp.csr_matrix(S)

                # Make A sparse
                A_sparse = sp.csr_matrix(A)

                # Set the b and beq vectors as numpy vectors
                b = np.array(b)
                beq = np.array(beq)

                # Add constraints
                model.addMConstr(Aeq_sparse, x, "=", beq, name="c")

                # Update the model to include the constraints added
                model.update()

                # Add constraints for the uneqalities of A
                model.addMConstr(A_sparse, x, "<", b, name="d")

                # Update the model with the extra constraints and then print it
                model.update()

                objective_function = np.asarray([-x for x in c])

                # Set the objective function in the model
                model.setMObjective(
                    None, objective_function, 0.0, None, None, x, GRB.MINIMIZE
                )
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
        print("Error code " + str(e.errno) + ": " + str(e))
    except AttributeError:
        print("Gurobi failed.")


def fast_fva(lb, ub, S, c, opt_percentage=100):
    """A Python function to perform fva using gurobi LP solver
    Returns the value of the optimal solution for all the following linear programs:
    min/max v_i, for all coordinates i=1,...,n, subject to,
    Sv = 0, lb <= v <= ub

    Keyword arguments:
    lb -- lower bounds for the fluxes, i.e., a n-dimensional vector
    ub -- upper bounds for the fluxes, i.e., a n-dimensional vector
    S -- the mxn stoichiometric matrix, s.t. Sv = 0
    c -- the objective function to maximize
    opt_percentage -- consider solutions that give you at least a certain
                      percentage of the optimal solution (default is to consider
                      optimal solutions only)
    """

    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )

    # declare the tolerance that gurobi works properly (we found it experimentally)
    tol = 1e-06

    m = S.shape[0]
    n = S.shape[1]
    beq = np.zeros(m)

    A = np.zeros((2 * n, n), dtype="float")
    A[0:n] = np.eye(n)
    A[n:] -= np.eye(n, n, dtype="float")

    b = np.concatenate((ub, -lb), axis=0)
    b = np.asarray(b, dtype="float")
    b = np.ascontiguousarray(b, dtype="float")

    # call fba to obtain an optimal solution
    max_biomass_flux_vector, max_biomass_objective = fast_fba(lb, ub, S, c)

    # add an additional constraint to impose solutions with at least `opt_percentage` of the optimal solution
    A = np.vstack((A, -c))

    b = np.append(
        b, -(opt_percentage / 100) * tol * math.floor(max_biomass_objective / tol)
    )

    min_fluxes = []
    max_fluxes = []

    try:

        # To avoid printint the output of the optimize() function of Gurobi, we need to set an environment like this
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()

            with gp.Model(env=env) as model:

                # Create variables
                x = model.addMVar(
                    shape=n,
                    vtype=GRB.CONTINUOUS,
                    name="x",
                    lb=-GRB.INFINITY,
                    ub=GRB.INFINITY,
                )

                # Make sparse Aeq
                Aeq_sparse = sp.csr_matrix(S)

                # Make A sparse
                A_sparse = sp.csr_matrix(A)

                # Set the b and beq vectors as numpy vectors
                b = np.array(b)
                beq = np.array(beq)

                # Add constraints
                model.addMConstr(Aeq_sparse, x, "=", beq, name="c")

                # Update the model to include the constraints added
                model.update()

                # Add constraints for the uneqalities of A
                model.addMConstr(A_sparse, x, "<", b, name="d")

                # Update the model with the extra constraints and then print it
                model.update()

                # Loop through the lines of the A matrix, set objective function for each and run the model
                for i in range(n):

                    # Set the ith row of the A matrix as the objective function
                    objective_function = A[
                        i,
                    ]

                    # Set the objective function in the model
                    model.setMObjective(
                        None, objective_function, 0.0, None, None, x, GRB.MINIMIZE
                    )
                    model.update()

                    # Optimize model
                    model.optimize()

                    # If optimized
                    status = model.status
                    if status == GRB.OPTIMAL:

                        # Get the min objective value
                        min_objective = model.getObjective().getValue()
                        min_fluxes.append(min_objective)
                    else:
                        min_fluxes.append(lb[i])

                    # Likewise, for the maximum
                    objective_function = np.asarray([-x for x in objective_function])
                    model.setMObjective(
                        None, objective_function, 0.0, None, None, x, GRB.MINIMIZE
                    )
                    model.update()
                    model.optimize()

                    # Again if optimized
                    status = model.status
                    if status == GRB.OPTIMAL:

                        # Get the max objective value
                        max_objective = -model.getObjective().getValue()
                        max_fluxes.append(max_objective)
                    else:
                        max_fluxes.append(ub[i])

                # Make lists of fluxes numpy arrays
                min_fluxes = np.asarray(min_fluxes)
                max_fluxes = np.asarray(max_fluxes)

                return (
                    min_fluxes,
                    max_fluxes,
                    max_biomass_flux_vector,
                    max_biomass_objective,
                )

    # Print error messages
    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))
    except AttributeError:
        print("Gurobi solver failed.")


def fast_find_redundant_facets(lb, ub, S, c, opt_percentage=100):
    """A Python function to perform fva using gurobi LP solver
    Returns the value of the optimal solution for all the following linear programs:
    min/max v_i, for all coordinates i=1,...,n, subject to,
    Sv = 0, lb <= v <= ub

    Keyword arguments:
    lb -- lower bounds for the fluxes, i.e., a n-dimensional vector
    ub -- upper bounds for the fluxes, i.e., a n-dimensional vector
    S -- the mxn stoichiometric matrix, s.t. Sv = 0
    c -- the objective function to maximize
    opt_percentage -- consider solutions that give you at least a certain
                      percentage of the optimal solution (default is to consider
                      optimal solutions only)
    """

    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )

    # declare the tolerance that gurobi works properly (we found it experimentally)
    red_facet_tol = 1e-07
    #shift_facet_tol = 1e-06
    tol = 1e-06

    m = S.shape[0]
    n = S.shape[1]
    beq = np.zeros(m)
    Aeq_res = S

    A = np.zeros((2 * n, n), dtype="float")
    A[0:n] = np.eye(n)
    A[n:] -= np.eye(n, n, dtype="float")

    b = np.concatenate((ub, -lb), axis=0)
    b = np.asarray(b, dtype="float")
    b = np.ascontiguousarray(b, dtype="float")

    # call fba to obtain an optimal solution
    max_biomass_flux_vector, max_biomass_objective = fast_fba(lb, ub, S, c)

    # add an additional constraint to impose solutions with at least `opt_percentage` of the optimal solution
    A = np.vstack((A, -c))

    b = np.append(
        b, -(opt_percentage / 100) * tol * math.floor(max_biomass_objective / tol)
    )

    b_res = []
    A_res = = np.zeros((0, 0), dtype="float")
    beq_res = np.array(beq)
    #Aeq_res = Aeq

    try:

        # To avoid printint the output of the optimize() function of Gurobi, we need to set an environment like this
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()

            with gp.Model(env=env) as model:

                # Create variables
                x = model.addMVar(
                    shape=n,
                    vtype=GRB.CONTINUOUS,
                    name="x",
                    lb=lb,
                    ub=ub,
                )

                # Make sparse Aeq
                Aeq_sparse = sp.csr_matrix(S)

                # Make A sparse
                A_sparse = sp.csr_matrix(A)

                # Set the b and beq vectors as numpy vectors
                b = np.array(b)
                beq = np.array(beq)

                # Add constraints
                model.addMConstr(Aeq_sparse, x, "=", beq, name="c")

                # Update the model to include the constraints added
                model.update()

                # Add constraints for the uneqalities of A
                #model.addMConstr(A_sparse, x, "<", b, name="d")

                # Update the model with the extra constraints and then print it
                #model.update()

                # Loop through the lines of the A matrix, set objective function for each and run the model
                for i in range(n):

                    # Set the ith row of the A matrix as the objective function
                    objective_function = A[
                        i,
                    ]

                    redundant_facet_right = False
                    redundant_facet_left = False

                    # for the maximum
                    objective_function_max = np.asarray([-x for x in objective_function])
                    model.setMObjective(
                        None, objective_function_max, 0.0, None, None, x, GRB.MINIMIZE
                    )
                    model.update()
                    model.optimize()

                    # Again if optimized
                    status = model.status
                    if status == GRB.OPTIMAL:
                        # Get the max objective value
                        max_objective = -model.getObjective().getValue()
                        #max_fluxes.append(max_objective)
                    else:
                        max_objective = ub[i]


                    model.variables[i].ub = ub[i]+10
                    model.update()
                    model.optimize()

                    status = model.status
                    if status == GRB.OPTIMAL:
                        # Get the max objective value
                        max_objective2 = -model.getObjective().getValue()
                        if np.abs(max_objective2 - max_objective) > red_facet_tol:
                            A_res = np.vstack(A_res, A[i,])
                            b_res.append(b[i])
                            redundant_facet_right = True




                    model.variables[i].ub = ub[i]
                    # for the minimum
                    #objective_function = np.asarray([-x for x in objective_function])
                    # Set the objective function in the model
                    model.setMObjective(
                        None, objective_function, 0.0, None, None, x, GRB.MINIMIZE
                    )
                    model.update()
                    # Optimize model
                    model.optimize()

                    # If optimized
                    status = model.status
                    if status == GRB.OPTIMAL:
                        # Get the min objective value
                        min_objective = model.getObjective().getValue()
                        #min_fluxes.append(min_objective)
                    else:
                        min_objective = lb[i]
                    

                    model.variables[i].lb = lb[i]-10
                    model.update()
                    model.optimize()

                    status = model.status
                    if status == GRB.OPTIMAL:
                        # Get the max objective value
                        min_objective2 = -model.getObjective().getValue()
                        if np.abs(min_objective2 - min_objective) > red_facet_tol:
                            A_res = np.vstack(A_res, A[n + i,])
                            b_res.append(b[n + i])
                            redundant_facet_left = True

                    model.variables[i].lb = lb[i]

                    if ((not redundant_facet_left) and (redundant_facet_right)):
                        width = abs(max_objective - min_objective)

                        # Check whether we keep or not the equality
                        if width < red_facet_tol:
                            Aeq_res = np.vstack(
                                (
                                    Aeq_res,
                                    A[
                                        i,
                                    ],
                                )
                            )
                            beq_res = np.append(beq_res, min(max_objective, min_objective))

                # Make lists of fluxes numpy arrays
                #min_fluxes = np.asarray(min_fluxes)
                #max_fluxes = np.asarray(max_fluxes)

                return (
                    A_res,
                    b_res,
                    Aeq_res,
                    beq_res,
                )

    # Print error messages
    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))
    except AttributeError:
        print("Gurobi solver failed.")



def fast_inner_ball(A, b):
    """A Python function to compute the maximum inscribed ball in the given polytope using gurobi LP solver
    Returns the optimal solution for the following linear program:
    max r, subject to,
    a_ix + r||a_i|| <= b, i=1,...,n

    Keyword arguments:
    A -- an mxn matrix that contains the normal vectors of the facets of the polytope row-wise
    b -- a m-dimensional vector
    """

    extra_column = []

    m = A.shape[0]
    n = A.shape[1]

    for i in range(A.shape[0]):
        entry = np.linalg.norm(
            A[
                i,
            ]
        )
        extra_column.append(entry)

    column = np.asarray(extra_column)
    A_expand = np.c_[A, column]

    with gp.Env(empty=True) as env:
        env.setParam("OutputFlag", 0)
        env.start()

        d = A_expand.shape[1]

        with gp.Model(env=env) as model:

            # Create variable
            x = model.addMVar(
                shape=d,
                vtype=GRB.CONTINUOUS,
                name="x",
                lb=-GRB.INFINITY,
                ub=GRB.INFINITY,
            )
            model.update()

            # Make A_full_dim sparse
            A_expand_sparse = sp.csr_matrix(A_expand.astype("float"))

            # Add constraints
            model.addMConstr(A_expand_sparse, x, "<", b, name="c")
            model.update()

            # Set the ith row of the A matrix as the objective function
            a = np.ones((1, n + 1))
            azero = np.zeros((1, n))
            a[:, :-1] = azero
            objective_function = a[0]

            # Set the objective function in the model
            model.setMObjective(
                None, objective_function, 0.0, None, None, x, GRB.MAXIMIZE
            )
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
            if r < 0:
                print(
                    "The radius calculated has negative value. The polytope is infeasible or something went wrong with the solver"
                )
            else:
                return point, r
