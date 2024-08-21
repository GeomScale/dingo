import pyoptinterface as poi
from pyoptinterface import highs, gurobi, copt, mosek
import numpy as np
import sys

default_solver = "highs"

def set_default_solver(solver_name):
    global default_solver
    default_solver = solver_name

def get_solver(solver_name):
    solvers = {"highs": highs, "gurobi": gurobi, "copt": copt, "mosek": mosek}
    if solver_name in solvers:
        return solvers[solver_name]
    else: 
        raise Exception("An unknown solver {solver_name} is requested.")

def dot(c, x):
    return poi.quicksum(c[i] * x[i] for i in range(len(x)) if abs(c[i]) > 1e-12)


def fba(lb, ub, S, c, solver_name=None):
    """A Python function to perform fba using PyOptInterface LP modeler
    Returns an optimal solution and its value for the following linear program:
    max c*v, subject to,
    Sv = 0, lb <= v <= ub

    Keyword arguments:
    lb -- lower bounds for the fluxes, i.e., a n-dimensional vector
    ub -- upper bounds for the fluxes, i.e., a n-dimensional vector
    S -- the mxn stoichiometric matrix, s.t. Sv = 0
    c -- the linear objective function, i.e., a n-dimensional vector
    solver_name  -- the solver to use, i.e., a string of the solver name (if None, the default solver is used)
    """

    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )
    if c.size != S.shape[1]:
        raise Exception(
            "The length of the linear objective function must be equal to the number of reactions."
        )

    m = S.shape[0]
    n = S.shape[1]
    optimum_value = 0
    optimum_sol = np.zeros(n)
    try:
        if solver_name is None:
            solver_name = default_solver
        SOLVER = get_solver(solver_name)
        # Create a model
        model = SOLVER.Model()
        model.set_model_attribute(poi.ModelAttribute.Silent, True)

        # Create variables and set lb <= v <= ub
        v = np.empty(n, dtype=object)
        for i in range(n):
            v[i] = model.add_variable(lb=lb[i], ub=ub[i])

        # Add the constraints Sv = 0
        for i in range(m):
            model.add_linear_constraint(dot(S[i], v), poi.Eq, 0)

        # Set the objective function
        obj = dot(c, v)
        model.set_objective(obj, sense=poi.ObjectiveSense.Maximize)

        # Optimize model
        model.optimize()

        # If optimized
        status = model.get_model_attribute(poi.ModelAttribute.TerminationStatus)
        if status == poi.TerminationStatusCode.OPTIMAL:
            optimum_value = model.get_value(obj)
            for i in range(n):
                optimum_sol[i] = model.get_value(v[i])
        return optimum_sol, optimum_value

    except poi.TerminationStatusCode.NUMERICAL_ERROR as e:
        print(f"A numerical error occurred: {e}")
    except poi.TerminationStatusCode.OTHER_ERROR as e:
        print(f"An error occurred: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def fva(lb, ub, S, c, opt_percentage=100, solver_name=None):
    """A Python function to perform fva using PyOptInterface LP modeler
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
    solver_name  -- the solver to use, i.e., a string of the solver name (if None, the default solver is used)
    """

    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )

    # declare the tolerance that highs and gurobi work properly (we found it experimentally)
    tol = 1e-06

    m = S.shape[0]
    n = S.shape[1]

    max_biomass_flux_vector, max_biomass_objective = fba(lb, ub, S, c, solver_name)

    min_fluxes = []
    max_fluxes = []

    adjusted_opt_threshold = (
        (opt_percentage / 100) * tol * np.floor(max_biomass_objective / tol)
    )

    try:
        if solver_name is None:
            solver_name = default_solver
        SOLVER = get_solver(solver_name)
        # Create a model
        model = SOLVER.Model()
        model.set_model_attribute(poi.ModelAttribute.Silent, True)

        # Create variables and set lb <= v <= ub
        v = np.empty(n, dtype=object)
        for i in range(n):
            v[i] = model.add_variable(lb=lb[i], ub=ub[i])

        # Add the constraints Sv = 0
        for i in range(m):
            model.add_linear_constraint(dot(S[i], v), poi.Eq, 0)

        # add an additional constraint to impose solutions with at least `opt_percentage` of the optimal solution
        model.add_linear_constraint(dot(c, v), poi.Geq, adjusted_opt_threshold)

        for i in range(n):
            # Set the objective function
            obj = poi.ExprBuilder(v[i])

            model.set_objective(obj, sense=poi.ObjectiveSense.Minimize)

            # Optimize model
            model.optimize()

            # If optimized
            status = model.get_model_attribute(poi.ModelAttribute.TerminationStatus)
            if status == poi.TerminationStatusCode.OPTIMAL:
                # Get the min objective value
                min_objective = model.get_value(v[i])
                min_fluxes.append(min_objective)
            else:
                min_fluxes.append(lb[i])

            # Likewise, for the maximum, optimize model
            model.set_objective(obj, sense=poi.ObjectiveSense.Maximize)

            # Again if optimized
            model.optimize()
            status = model.get_model_attribute(poi.ModelAttribute.TerminationStatus)
            if status == poi.TerminationStatusCode.OPTIMAL:
                # Get the max objective value
                max_objective = model.get_value(v[i])
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

    except poi.TerminationStatusCode.NUMERICAL_ERROR as e:
        print(f"A numerical error occurred: {e}")
    except poi.TerminationStatusCode.OTHER_ERROR as e:
        print(f"An error occurred: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def inner_ball(A, b, solver_name=None):
    """A Python function to compute the maximum inscribed ball in the given polytope using PyOptInterface LP modeler
    Returns the optimal solution for the following linear program:
    max r, subject to,
    a_ix + r||a_i|| <= b, i=1,...,n

    Keyword arguments:
    A -- an mxn matrix that contains the normal vectors of the facets of the polytope row-wise
    b -- a m-dimensional vector
    solver_name  -- the solver to use, i.e., a string of the solver name (if None, the default solver is used)
    """

    extra_column = []

    m = A.shape[0]
    n = A.shape[1]

    for i in range(A.shape[0]):
        entry = np.linalg.norm(A[i])
        extra_column.append(entry)

    column = np.asarray(extra_column)
    A_expand = np.c_[A, column]

    if solver_name is None:
        solver_name = default_solver
    SOLVER = get_solver(solver_name)
    model = SOLVER.Model()
    model.set_model_attribute(poi.ModelAttribute.Silent, True)

    # Create variables where x[n] is the radius
    x = np.empty(n + 1, dtype=object)
    for i in range(n + 1):
        x[i] = model.add_variable()

    # Add the constraints a_ix + r||a_i|| <= b
    for i in range(m):
        model.add_linear_constraint(dot(A_expand[i], x), poi.Leq, b[i])

    # Set the objective function
    obj = poi.ExprBuilder(x[n])
    model.set_objective(obj, sense=poi.ObjectiveSense.Maximize)

    # Optimize model
    model.optimize()

    # Get the center point and the radius of max ball from the solution of LP
    point = [model.get_value(x[i]) for i in range(n)]

    # Get radius
    r = model.get_value(obj)

    # And check whether the computed radius is negative
    if r < 0:
        raise Exception(
            "The radius calculated has negative value. The polytope is infeasible or something went wrong with the solver"
        )
    else:
        return point, r


def set_model(n, lb, ub, Aeq, beq, A, b, solver_name=None):
    """
    A helper function of remove_redundant_facets function
    Create a PyOptInterface model with given PyOptInterface variables, equality constraints, inequality constraints and solver name
    but without an objective function.
    """
    # Create a model
    if solver_name is None:
        solver_name = default_solver
    SOLVER = get_solver(solver_name)
    model = SOLVER.Model()
    model.set_model_attribute(poi.ModelAttribute.Silent, True)

    # Create variables
    x = np.empty(n, dtype=object)
    for i in range(n):
        x[i] = model.add_variable(lb=lb[i], ub=ub[i])

    # Add the equality constraints
    for i in range(Aeq.shape[0]):
        model.add_linear_constraint(dot(Aeq[i], x), poi.Eq, beq[i])

    # Add the inequality constraints
    for i in range(A.shape[0]):
        model.add_linear_constraint(dot(A[i], x), poi.Leq, b[i])

    return model, x


def remove_redundant_facets(lb, ub, S, c, opt_percentage=100, solver_name=None):
    """A function to find and remove the redundant facets and to find
    the facets with very small offset and to set them as equalities

    Keyword arguments:
    lb -- lower bounds for the fluxes, i.e., a n-dimensional vector
    ub -- upper bounds for the fluxes, i.e., a n-dimensional vector
    S -- the mxn stoichiometric matrix, s.t. Sv = 0
    c -- the objective function to maximize
    opt_percentage -- consider solutions that give you at least a certain
                      percentage of the optimal solution (default is to consider
                      optimal solutions only)
    solver_name  -- the solver to use, i.e., a string of the solver name (if None, the default solver is used)
    """

    if lb.size != S.shape[1] or ub.size != S.shape[1]:
        raise Exception(
            "The number of reactions must be equal to the number of given flux bounds."
        )

    # declare the tolerance that highs and gurobi work properly (we found it experimentally)
    redundant_facet_tol = 1e-07
    tol = 1e-06

    m = S.shape[0]
    n = S.shape[1]

    # [v,-v] <= [ub,-lb]
    A = np.zeros((2 * n, n), dtype="float")
    A[0:n] = np.eye(n)
    A[n:] -= np.eye(n, n, dtype="float")

    b = np.concatenate((ub, -lb), axis=0)
    b = np.ascontiguousarray(b, dtype="float")

    beq = np.zeros(m)

    Aeq_res = S
    beq_res = np.array(beq)
    b_res = []
    A_res = np.empty((0, n), float)

    max_biomass_flux_vector, max_biomass_objective = fba(lb, ub, S, c, solver_name)
    val = -np.floor(max_biomass_objective / tol) * tol * opt_percentage / 100

    start = np.zeros(n)

    try:

        # initialize
        indices_iter = range(n)
        removed = 1
        offset = 1
        facet_left_removed = np.zeros(n, dtype=bool)
        facet_right_removed = np.zeros(n, dtype=bool)

        # Loop until no redundant facets are found
        while removed > 0 or offset > 0:
            removed = 0
            offset = 0
            indices = indices_iter
            indices_iter = []

            Aeq = np.array(Aeq_res)
            beq = np.array(beq_res)

            A_res = np.empty((0, n), dtype=float)
            b_res = []

            model_iter, v = set_model(n, lb, ub, Aeq, beq, np.array([-c]), [val], solver_name)

            for cnt, i in enumerate(indices):

                redundant_facet_right = True
                redundant_facet_left = True

                if cnt > 0:
                    last_idx = indices[cnt-1]
                    model_iter.set_variable_attribute(
                        v[last_idx], poi.VariableAttribute.LowerBound, lb[last_idx]
                    )
                    model_iter.set_variable_attribute(
                        v[last_idx], poi.VariableAttribute.UpperBound, ub[last_idx]
                    )
                    
                # objective function
                obj = poi.ExprBuilder(v[i])

                # maximize v_i (right)
                model_iter.set_objective(obj, sense=poi.ObjectiveSense.Maximize)
                model_iter.optimize()

                # if optimized
                status = model_iter.get_model_attribute(
                    poi.ModelAttribute.TerminationStatus
                )
                if status == poi.TerminationStatusCode.OPTIMAL:
                    # get the maximum objective value
                    max_objective = model_iter.get_value(obj)
                else:
                    max_objective = ub[i]

                # if this facet was not removed in a previous iteration
                if not facet_right_removed[i]:
                    # Relax the inequality
                    model_iter.set_variable_attribute(
                        v[i], poi.VariableAttribute.UpperBound, ub[i] + 1
                    )

                    # Solve the model
                    model_iter.optimize()

                    status = model_iter.get_model_attribute(
                        poi.ModelAttribute.TerminationStatus
                    )
                    if status == poi.TerminationStatusCode.OPTIMAL:
                        # Get the max objective value with relaxed inequality

                        max_objective2 = model_iter.get_value(obj)
                        if np.abs(max_objective2 - max_objective) > redundant_facet_tol:
                            redundant_facet_right = False
                        else:
                            removed += 1
                            facet_right_removed[i] = True

                    # Reset the inequality
                    model_iter.set_variable_attribute(
                        v[i], poi.VariableAttribute.UpperBound, ub[i]
                    )

                # minimum v_i (left)
                model_iter.set_objective(obj, sense=poi.ObjectiveSense.Minimize)
                model_iter.optimize()

                # If optimized
                status = model_iter.get_model_attribute(
                    poi.ModelAttribute.TerminationStatus
                )
                if status == poi.TerminationStatusCode.OPTIMAL:
                    # Get the min objective value
                    min_objective = model_iter.get_value(obj)
                else:
                    min_objective = lb[i]

                # if this facet was not removed in a previous iteration
                if not facet_left_removed[i]:
                    # Relax the inequality
                    model_iter.set_variable_attribute(
                        v[i], poi.VariableAttribute.LowerBound, lb[i] - 1
                    )

                    # Solve the model
                    model_iter.optimize()

                    status = model_iter.get_model_attribute(
                        poi.ModelAttribute.TerminationStatus
                    )
                    if status == poi.TerminationStatusCode.OPTIMAL:
                        # Get the min objective value with relaxed inequality
                        min_objective2 = model_iter.get_value(obj)
                        if np.abs(min_objective2 - min_objective) > redundant_facet_tol:
                            redundant_facet_left = False
                        else:
                            removed += 1
                            facet_left_removed[i] = True

                if (not redundant_facet_left) or (not redundant_facet_right):
                    width = abs(max_objective - min_objective)

                    # Check whether the offset in this dimension is small (and set an equality)
                    if width < redundant_facet_tol:
                        offset += 1
                        Aeq_res = np.vstack((Aeq_res, A[i]))
                        beq_res = np.append(beq_res, min(max_objective, min_objective))
                        # Remove the bounds on this dimension
                        ub[i] = sys.float_info.max
                        lb[i] = -sys.float_info.max
                    else:
                        # store this dimension
                        indices_iter.append(i)

                        if not redundant_facet_left:
                            # Not a redundant inequality
                            A_res = np.append(A_res, np.array([A[n + i]]), axis=0)
                            b_res.append(b[n + i])
                        else:
                            lb[i] = -sys.float_info.max

                        if not redundant_facet_right:
                            # Not a redundant inequality
                            A_res = np.append(A_res, np.array([A[i]]), axis=0)
                            b_res.append(b[i])
                        else:
                            ub[i] = sys.float_info.max
                else:
                    # Remove the bounds on this dimension
                    ub[i] = sys.float_info.max
                    lb[i] = -sys.float_info.max

        b_res = np.asarray(b_res, dtype="float")
        A_res = np.asarray(A_res, dtype="float")
        A_res = np.ascontiguousarray(A_res, dtype="float")
        return A_res, b_res, Aeq_res, beq_res

    except poi.TerminationStatusCode.NUMERICAL_ERROR as e:
        print(f"A numerical error occurred: {e}")
    except poi.TerminationStatusCode.OTHER_ERROR as e:
        print(f"An error occurred: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
