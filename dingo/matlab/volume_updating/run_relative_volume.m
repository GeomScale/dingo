function result = run_relative_volume(model_path, lambda, user_opts)
    %RUN_RELATIVE_VOLUME Load a metabolic model, generate a random halfspace,
    % and compute the relative volume vol(P ∩ H(c0)) / vol(P).
    %
    % Input:
    %   model_path - path to a .mat file containing a 'problem' struct
    %                with fields .Aeq (=S), .beq, .lb, .ub.
    %                Relative to repo root, e.g. 'data/e_coli_core'
    %                (the .mat extension is optional).
    %   lambda     - scalar in (0,1) determining the halfspace threshold:
    %                  c0 = (1-lambda)*c_max + lambda*c_min
    %                lambda=0  → c0 = c_max (halfspace contains all of P)
    %                lambda=1  → c0 = c_min (halfspace barely touches P)
    %   user_opts  - optional struct overriding default parameters for
    %                relative_volume (see default_params.m).
    %
    % Output:
    %   result - struct containing all fields from relative_volume plus:
    %     .rho      the relative volume estimate
    %     .c        random cost vector used (d x 1, unit norm)
    %     .c_min    min c'*v over P (from LP)
    %     .c_max    max c'*v over P (from LP)
    %     .c0       halfspace threshold c'*v <= c0
    %     .lambda   the lambda value used

    % -----------------------------------------------------------------
    % 0. Path setup & load the model
    % -----------------------------------------------------------------
    avail_simdLen = init();  %#ok<NASGU>  % sets up paths, compiles solvers

    root = fileparts(fileparts(mfilename('fullpath')));

    % Accept path with or without .mat extension
    if ~endsWith(model_path, '.mat')
        full_path = fullfile(root, [model_path '.mat']);
    else
        full_path = fullfile(root, model_path);
    end

    if ~exist(full_path, 'file')
        error('Model file not found: %s', full_path);
    end

    data = load(full_path);

    if isfield(data, 'problem')
        % --- Existing format: file contains a 'problem' struct ---
        problem = data.problem;
    else
        % --- Try to auto-detect a COBRA-style metabolic model ---
        model = detect_metabolic_model(data);
        if isempty(model)
            error('File does not contain a ''problem'' struct or a recognized metabolic model: %s', full_path);
        end
        fprintf('  Detected metabolic model: %d rxns, %d mets\n', ...
                length(model.rxns), length(model.mets));
        problem = metabolic_to_problem(model);
    end
    problem = standardize_problem(problem);

    % -----------------------------------------------------------------
    % Detect and handle netlib-style (unbounded) problems.
    % Netlib .mat files ship with a non-empty df vector (the LP cost).
    % Metabolic models have df = [] (uniform).  Netlib polytopes are
    % typically unbounded, so we apply the same bounding logic used by
    % PolytopeSampler's loadProblem.m.
    % -----------------------------------------------------------------
    is_netlib = ~isempty(problem.df) && isnumeric(problem.df);
    if is_netlib
        fprintf('  Detected netlib-style problem (non-empty df).  Applying bounds...\n');

        Aeq_net = problem.Aeq;
        beq_net = problem.beq(:);
        lb_net  = problem.lb(:);
        ub_net  = problem.ub(:);
        if isfield(problem, 'Aineq') && ~isempty(problem.Aineq)
            Aineq_net = problem.Aineq;
            bineq_net = problem.bineq(:);
        else
            Aineq_net = [];
            bineq_net = [];
        end

        lp_opts = optimoptions('linprog', 'Display', 'none');

        % Solve the netlib LP:  min  df'*x  over the original constraints
        df_net = problem.df(:);
        [x_opt, ~, exitflag] = linprog(df_net, Aineq_net, bineq_net, ...
                                        Aeq_net, beq_net, lb_net, ub_net, lp_opts);

        if exitflag <= 0
            warning('Netlib bounding LP did not converge (exitflag=%d). Using problem as-is.', exitflag);
        else
            % Add a bounding inequality:  df'*x <= df'*x_opt + |df|'*|x_opt|
            threshold = df_net' * x_opt + abs(df_net)' * abs(x_opt);
            problem.Aineq = [Aineq_net; df_net'];
            problem.bineq = [bineq_net; threshold];

            % Clamp variable bounds around the optimal vertex
            bound = 2 * max(abs(x_opt));
            problem.lb = max(lb_net, -bound);
            problem.ub = min(ub_net,  bound);

            % Switch to uniform sampling
            problem.f   = [];
            problem.df  = [];
            problem.ddf = [];

            fprintf('    Bounding LP solved.  Added df''*x <= %.4f, bounds clamped to ±%.4f.\n', ...
                    threshold, bound);
        end
    end

    % Effective dimension (before presolve — number of original variables)
    d = max([size(problem.Aeq, 2), length(problem.lb), length(problem.ub)]);

    % -----------------------------------------------------------------
    % 1. Generate random cost vector
    % -----------------------------------------------------------------
    c = randn(d, 1);
    c = c / norm(c);

    % -----------------------------------------------------------------
    % 2. Compute c_min and c_max via LP
    % -----------------------------------------------------------------
    Aeq = problem.Aeq;
    beq = problem.beq(:);
    lb  = problem.lb(:);
    ub  = problem.ub(:);

    if isfield(problem, 'Aineq') && ~isempty(problem.Aineq)
        Aineq = problem.Aineq;
        bineq = problem.bineq(:);
    else
        Aineq = [];
        bineq = [];
    end

    lp_opts = optimoptions('linprog', 'Display', 'none');

    % min  c'*x  over P
    [~, c_min] = linprog(c, Aineq, bineq, Aeq, beq, lb, ub, lp_opts);

    % max  c'*x = -min -c'*x  over P
    [~, c_max_neg] = linprog(-c, Aineq, bineq, Aeq, beq, lb, ub, lp_opts);
    c_max = -c_max_neg;

    % -----------------------------------------------------------------
    % 3. Set threshold: interpolate between c_min and c_max
    % -----------------------------------------------------------------
    c0 = (1 - lambda) * c_max + lambda * c_min;

    fprintf('=== run_relative_volume ===\n');
    fprintf('  model     : %s\n', model_path);
    fprintf('  dimension : %d\n', d);
    fprintf('  lambda    : %.4f\n', lambda);
    fprintf('  c_min     : %.4f\n', c_min);
    fprintf('  c_max     : %.4f\n', c_max);
    fprintf('  c0        : %.4f\n', c0);
    fprintf('  (lambda=0 => c0=c_max=%.4f, lambda=1 => c0=c_min=%.4f)\n', c_max, c_min);

    % -----------------------------------------------------------------
    % 4. Call relative_volume
    % -----------------------------------------------------------------
    if nargin < 3 || isempty(user_opts)
        [rho, vol_result] = relative_volume(problem, c, c0);
    else
        [rho, vol_result] = relative_volume(problem, c, c0, user_opts);
    end

    % -----------------------------------------------------------------
    % 5. Assemble output
    % -----------------------------------------------------------------
    result = vol_result;
    result.rho    = rho;
    result.c      = c;
    result.c_min  = c_min;
    result.c_max  = c_max;
    result.c0     = c0;
    result.lambda = lambda;
end

% =========================================================================
%  Helper: detect a COBRA-style metabolic model in a loaded .mat struct
% =========================================================================
function model = detect_metabolic_model(data)
    %DETECT_METABOLIC_MODEL Find the first COBRA-style model in loaded data.
    %
    % A metabolic model is identified by having fields: S, lb, ub, and rxns.

    fields = fieldnames(data);
    for i = 1:numel(fields)
        var = data.(fields{i});
        if isstruct(var) && all(isfield(var, {'S', 'lb', 'ub', 'rxns'}))
            % Found a COBRA model — store its variable name for reporting
            model = var;
            model.model_var_name = fields{i};  %#ok<AGROW>
            return;
        end
    end

    % No metabolic model found
    model = [];
end

% =========================================================================
%  Helper: convert a COBRA metabolic model to the PolytopeSampler problem
%          struct (Aeq·x = beq, lb ≤ x ≤ ub).
% =========================================================================
function problem = metabolic_to_problem(model)
    %METABOLIC_TO_PROBLEM Convert COBRA model to PolytopeSampler format.
    %
    % COBRA model fields used:
    %   .S   - stoichiometric matrix (mets × rxns),  S·v = b
    %   .b   - right-hand side (default: zeros)
    %   .lb  - lower bounds on reaction fluxes
    %   .ub  - upper bounds on reaction fluxes
    %
    % Output problem fields:
    %   .Aeq = S  (after removing all-zero metabolite rows)
    %   .beq = b  (after removing the same rows)
    %   .lb  = lb(:)
    %   .ub  = ub(:)
    %   .f, .df, .ddf = []   (uniform sampling)

    % --- Stoichiometric matrix & right-hand side ---
    S = model.S;
    if isfield(model, 'b') && ~isempty(model.b)
        b_vec = model.b(:);
    else
        b_vec = zeros(size(S, 1), 1);
    end

    % Remove metabolite rows that impose no constraint (all-zero rows)
    nonzero_rows = any(S, 2);
    S = S(nonzero_rows, :);
    b_vec = b_vec(nonzero_rows);

    % --- Assemble problem struct ---
    problem.Aeq = S;
    problem.beq = b_vec;
    problem.lb  = model.lb(:);
    problem.ub  = model.ub(:);
    problem.f   = [];
    problem.df  = [];
    problem.ddf = [];
end
