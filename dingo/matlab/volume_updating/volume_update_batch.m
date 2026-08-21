function volume_update_batch(request_path, output_path)
%VOLUME_UPDATE_BATCH Apply sequential cuts using relative_volume unchanged.
%
% The request MAT file must contain S, lb, ub, cut_normals and
% cut_thresholds.  Each new cut becomes part of the parent polytope before the
% next call.  No objective, biomass floor, FBA or alternative volume logic is
% introduced here.

    req = load(request_path);
    required = {'S', 'lb', 'ub', 'cut_normals', 'cut_thresholds'};
    for k = 1:numel(required)
        if ~isfield(req, required{k})
            error('Missing request field: %s', required{k});
        end
    end

    S = double(req.S);
    lb = double(req.lb(:));
    ub = double(req.ub(:));
    cut_normals = double(req.cut_normals);
    cut_thresholds = double(req.cut_thresholds(:));
    n = numel(lb);

    if size(S, 2) ~= n || numel(ub) ~= n
        error('Dimension mismatch among S, lb and ub.');
    end
    if size(cut_normals, 2) ~= n
        error('cut_normals must have one column per reaction.');
    end
    if size(cut_normals, 1) ~= numel(cut_thresholds)
        error('cut_normals rows and cut_thresholds length differ.');
    end
    if any(~isfinite(S(:))) || any(~isfinite(lb)) || any(~isfinite(ub)) || ...
       any(~isfinite(cut_normals(:))) || any(~isfinite(cut_thresholds))
        error('All wrapper inputs must be finite.');
    end
    if any(lb > ub)
        error('At least one lower bound exceeds its upper bound.');
    end

    problem.Aeq = S;
    problem.beq = zeros(size(S, 1), 1);
    problem.Aineq = zeros(0, n);
    problem.bineq = zeros(0, 1);
    problem.lb = lb;
    problem.ub = ub;
    problem.f = [];
    problem.df = [];
    problem.ddf = [];

    if isfield(req, 'options')
        options = req.options;
    else
        options = struct();
    end

    n_cuts = size(cut_normals, 1);
    cut_ratios = ones(n_cuts, 1);
    cut_log_ratios = zeros(n_cuts, 1);
    cut_n_phases = zeros(n_cuts, 1);
    cut_total_steps = zeros(n_cuts, 1);
    cut_converged = false(n_cuts, 1);
    cut_diagnostics = cell(n_cuts, 1);

    for j = 1:n_cuts
        c = cut_normals(j, :)';
        c0 = cut_thresholds(j);

        [rho_j, diag_j] = relative_volume(problem, c, c0, options);
        if ~isscalar(rho_j) || ~isfinite(rho_j) || rho_j < 0 || rho_j > 1 + 1e-6
            error('Invalid ratio returned for cut %d: %.17g', j, rho_j);
        end

        cut_ratios(j) = rho_j;
        cut_log_ratios(j) = log(max(rho_j, 1e-300));
        if isfield(diag_j, 'n_phases'); cut_n_phases(j) = diag_j.n_phases; end
        if isfield(diag_j, 'total_steps'); cut_total_steps(j) = diag_j.total_steps; end
        if isfield(diag_j, 'converged')
            cut_converged(j) = all(diag_j.converged(:));
        end
        % MATLAB handle objects cannot be serialized reliably to a portable
        % v7 MAT file.  They are internal sampler state, not wrapper output.
        if isfield(diag_j, 'samplers')
            diag_j = rmfield(diag_j, 'samplers');
        end
        cut_diagnostics{j} = diag_j;

        % This is the entire sequential-update rule: the accepted cut becomes
        % part of the next parent, and the same static Apostolos routine is
        % called again for the next cut.
        problem.Aineq(end + 1, :) = c'; %#ok<AGROW>
        problem.bineq(end + 1, 1) = c0; %#ok<AGROW>
    end

    log_rho = sum(cut_log_ratios);
    rho = exp(log_rho);
    algorithm = 'static_relative_volume_sequential';
    objective_floor = false;
    save(output_path, 'rho', 'log_rho', 'cut_ratios', 'cut_log_ratios', ...
         'cut_n_phases', 'cut_total_steps', 'cut_converged', ...
         'cut_diagnostics', 'algorithm', 'objective_floor', '-v7');
end
