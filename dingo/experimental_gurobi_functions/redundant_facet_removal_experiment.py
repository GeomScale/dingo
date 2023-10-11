# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis

# Licensed under GNU LGPL.3, see LICENCE file

import sys
import numpy as np
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB
import math

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


def solve_lp_with_different_objectives(model, new_objective_coeffs):
    """
    Solve a linear program with a different objective function.

    Parameters:
        model (gurobipy.Model): The Gurobi model with the original constraints.
        new_objective_coeffs (list or numpy.ndarray): List or array of new objective coefficients for the variables.

    Returns:
        gurobipy.Model: The updated Gurobi model with the new objective function.
    """
    # Clear the existing objective function
    model.setObjective(0, clear=True)

    # TODO: Optimization opportunity - You could pass a numpy array here (instead of a list of coefficients)
    # and update the objective in one call without a for loop.
    # For example, you can use:
    # model.setObjective(new_objective_coeffs, GRB.MINIMIZE)

    # Update the objective function with the new coefficients
    for i, var in enumerate(model.getVars()):
        var.setAttr(GRB.Attr.Obj, new_objective_coeffs[i])

    # Optimize the updated model
    model.optimize()

    return model


def fast_remove_redundant_facets(lb, ub, S, c, opt_percentage=100):
    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )

    redundant_facet_tol = 1e-07
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

    max_biomass_flux_vector, max_biomass_objective = fast_fba(lb, ub, S, c)
    val = -np.floor(max_biomass_objective / tol) * tol * opt_percentage / 100

    b_res = []
    A_res = np.empty((0, n), float)
    beq_res = np.array(beq)

    try:
        with gp.Env(empty=True) as env:
            env.setParam("OutputFlag", 0)
            env.start()

            with gp.Model(env=env) as model:
                x = model.addMVar(
                    shape=n,
                    vtype=GRB.CONTINUOUS,
                    name="x",
                    lb=lb,
                    ub=ub,
                )

                Aeq_sparse = sp.csr_matrix(S)
                A_sparse = sp.csr_matrix(np.array(-c))
                b_sparse = np.array(val)

                b = np.array(b)
                beq = np.array(beq)

                model.addMConstr(Aeq_sparse, x, "=", beq, name="c")
                model.update()
                model.addMConstr(A_sparse, x, "<", [val], name="d")
                model.update()

                model_iter = model.copy()

                indices_iter = range(n)
                removed = 1
                offset = 1
                facet_left_removed = np.zeros((1, n), dtype=bool)
                facet_right_removed = np.zeros((1, n), dtype=bool)

                while removed > 0 or offset > 0:
                    removed = 0
                    offset = 0
                    indices = indices_iter
                    indices_iter = []

                    Aeq_sparse = sp.csr_matrix(Aeq_res)
                    beq = np.array(beq_res)

                    b_res = []
                    A_res = np.empty((0, n), float)
                   
                    update_model_constraints_and_bounds(model_iter, Aeq_sparse, beq, A_sparse, [val], lb, ub)


                    for i in indices:
                        objective_function = A[i]

                        redundant_facet_right = True
                        redundant_facet_left = True

                        objective_function_max = np.asarray(
                            [-x for x in objective_function]
                        )

                        model_iter = solve_lp_with_different_objectives(
                            model_iter.copy(), objective_function_max
                        )


                        status = model_iter.status
                        if status == GRB.OPTIMAL:
                            max_objective = -model_iter.objVal
                        else:
                            max_objective = ub[i]

                        if not facet_right_removed[0, i]:
                            ub_iter = ub.copy()
                            ub_iter[i] = ub_iter[i] + 1

                            # Call solve_lp_with_different_objectives to solve LP
                            model_iter = solve_lp_with_different_objectives(
                                model_iter.copy(), objective_function
                            )

                            status = model_iter.status
                            if status == GRB.OPTIMAL:
                                max_objective2 = -model_iter.objVal
                                if (
                                    np.abs(max_objective2 - max_objective)
                                    > redundant_facet_tol
                                ):
                                    redundant_facet_right = False
                                else:
                                    removed += 1
                                    facet_right_removed[0, i] = True

                        model_iter.reset()
                        x = model_iter.getVars()
                        for j in range(n):
                            x[j].LB = lb[j]
                            x[j].UB = ub[j]
                            x[j].obj = objective_function[j]

                        model_iter.optimize()

                        status = model_iter.status
                        if status == GRB.OPTIMAL:
                            min_objective = model_iter.objVal
                        else:
                            min_objective = lb[i]

                        if not facet_left_removed[0, i]:
                            lb_iter = lb.copy()
                            lb_iter[i] = lb_iter[i] - 1

                            # Call solve_lp_with_different_objectives to solve LP
                            model_iter = solve_lp_with_different_objectives(
                                model_iter.copy(), objective_function
                            )


                            status = model_iter.status
                            if status == GRB.OPTIMAL:
                                min_objective2 = model_iter.objVal
                                if (
                                    np.abs(min_objective2 - min_objective)
                                    > redundant_facet_tol
                                ):
                                    redundant_facet_left = False
                                else:
                                    removed += 1
                                    facet_left_removed[0, i] = True

                        if (not redundant_facet_left) or (not redundant_facet_right):
                            width = abs(max_objective - min_objective)

                            if width < redundant_facet_tol:
                                offset += 1
                                Aeq_res = np.vstack(
                                    (
                                        Aeq_res,
                                        A[
                                            i,
                                        ],
                                    )
                                )
                                beq_res = np.append(
                                    beq_res, min(max_objective, min_objective)
                                )
                                ub[i] = sys.float_info.max
                                lb[i] = -sys.float_info.max
                            else:
                                indices_iter.append(i)

                                if not redundant_facet_left:
                                    A_res = np.append(
                                        A_res,
                                        np.array(
                                            [
                                                A[
                                                    n + i,
                                                ]
                                            ]
                                        ),
                                        axis=0,
                                    )
                                    b_res.append(b[n + i])
                                else:
                                    lb[i] = -sys.float_info.max

                                if not redundant_facet_right:
                                    A_res = np.append(
                                        A_res,
                                        np.array(
                                            [
                                                A[
                                                    i,
                                                ]
                                            ]
                                        ),
                                        axis=0,
                                    )
                                    b_res.append(b[i])
                                else:
                                    ub[i] = sys.float_info.max
                        else:
                            ub[i] = sys.float_info.max
                            lb[i] = -sys.float_info.max

                b_res = np.asarray(b_res)
                A_res = np.asarray(A_res, dtype="float")
                A_res = np.ascontiguousarray(A_res, dtype="float")

                return A_res, b_res, Aeq_res, beq_res

    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))
    except AttributeError:
        print("Gurobi solver failed.")


