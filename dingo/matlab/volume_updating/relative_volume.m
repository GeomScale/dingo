function [rho, result] = relative_volume(problem, c, c0, user_opts)
    %RELATIVE_VOLUME Estimate vol(P ∩ H(c0)) / vol(P) via halfspace simulated annealing.
    %
    % Input:
    %   problem - struct with fields .Aeq (=S), .beq (=0), .lb, .ub
    %             (and optionally .Aineq, .bineq) defining the metabolic polytope
    %             P = {v | Sv=0, lb <= v <= ub}.
    %   c       - cost vector (d x 1 or 1 x d)
    %   c0      - scalar target threshold for halfspace c'*v <= c0
    %   user_opts - optional struct overriding default parameters (see default_params.m)
    %
    % Output:
    %   rho    - estimate of vol(P ∩ H(c0)) / vol(P)
    %   result - struct with diagnostic fields:
    %     .thresholds     cell array of thresholds [t_max_or_Inf; t_k; ...; t_1]
    %     .ratios         individual ratio estimates [r_k, ..., r_0]
    %     .log_ratios     log of each ratio estimate
    %     .steps_per_phase  HMC step() counts per ratio
    %     .total_steps    total HMC step() calls
    %     .n_phases       number of intermediate bodies (= length(thresholds)-1)
    %     .log_rho        log of final estimate
    %     .converged      whether sliding window converged for each ratio
    %
    % Algorithm: Chalkis et al. (2023) annealing schedule with
    %   Cousins & Vempala (2016) sliding window convergence,
    %   using Kook et al. (2022) CRHMC via PolytopeSamplerMatlab.

    % =====================================================================
    % 0. Input validation and setup
    % =====================================================================
    c = c(:);  % ensure column vector

    % Shuffle the global RNG so each run produces different seeds / chains
    rng('shuffle');

    % Must initialize paths BEFORE using Setfield or default_options.
    % init() returns the maximum usable simdLen (8 if MEX compiled, 4 if
    % only pre-compiled binaries are available).
    avail_simdLen = init();

    if nargin < 4 || isempty(user_opts)
        params = default_params();
    else
        params = Setfield(default_params(), user_opts);
    end
    params

    % Use the best simdLen that is actually available
    if params.simdLen > avail_simdLen
        if params.verb >= 1
            fprintf('Note: simdLen reduced from %d to %d (MEX solver unavailable).\n', ...
                    params.simdLen, avail_simdLen);
        end
        params.simdLen = avail_simdLen;
    end

    [polytope_opts, sampler_opts] = build_opts(params);

    % =====================================================================
    % Edge case: empty polytope or invalid problem
    % =====================================================================
    n_vars = max([size(problem.Aeq,2), length(problem.lb), length(problem.ub)]);
    if n_vars == 0
        error('Problem has zero variables.');
    end

    % =====================================================================
    % Phase 1: Build the annealing schedule
    % =====================================================================
    if params.verb >= 1
        fprintf('=== Phase 1: Building Annealing Schedule ===\n');
    end

    % Create full polytope P (no extra halfspace)
    registry = struct();
    registry.thresholds = {};
    registry.polytopes  = {};
    registry.samplers   = {};
    registry.initial_ratios = {};   % warm-start values for Phase 2
    registry.dim        = 0;
    registry.c          = c;

    registry = registry_store(registry, Inf, c, problem, polytope_opts, sampler_opts);
    sampler_full = registry.samplers{1};
    registry.dim = size(sampler_full.x, 2);

    % ---- Phase 1A+B: Collect first batch, determine t_max, find t_star ----
    % Draw nu*N points from the full polytope P.  These serve double duty:
    %   (a) determine t_max = max_i (c'*x_i) * safety
    %   (b) reuse as the first U/L-test batch in the binary search for t_star
    if params.verb >= 1
        fprintf('  Sampling %d points from P to estimate t_max...\n', ...
                params.nu * params.N_utest);
    end

    n_needed = params.nu * params.N_utest;
    [points, sampler_full] = advance_sampler(sampler_full, n_needed);
    c_vals = points * c;
    t_max = max(c_vals) * params.t_max_safety;

    if params.verb >= 1
        fprintf('    t_max = %.4f (safety factor %.2f)\n', t_max, params.t_max_safety);
        fprintf('    target c0 = %.4f\n', c0);
    end

    % Edge case: halfspace contains all of P
    if c0 >= t_max
        if params.verb >= 1
            fprintf('  c0 >= t_max: halfspace contains entire polytope. rho = 1.\n');
        end
        rho = 1.0;
        result = build_result(registry, [], [], [], 0, 0, true(0,1));
        return;
    end

    % ---- Check: is the ratio to P0 already large enough? ----
    % Run the L-test on the initial points: H0: theta <= r  vs  H1: theta > r.
    % If it passes, vol(P0)/vol(P) >= r with confidence 1-alpha, i.e. the
    % ratio is already large enough for direct rejection sampling.  Skip the
    % binary search and all intermediate body construction.
    [~, L_pass_direct, theta_bar_direct, ~] = ul_test(points, c0, c, params);

    if params.verb >= 1
        fprintf('    L-test to P0 on initial points: theta_bar=%.4f, L_pass=%d\n', ...
                theta_bar_direct, L_pass_direct);
    end

    if L_pass_direct
        if params.verb >= 1
            fprintf('    Ratio is already >= r — direct estimation (1 phase).\n');
        end
        % sampler_full is the only sampler we need; Phase 2 estimates
        % vol(P0) / vol(P) directly via rejection sampling.
        direct_estimation = true;
        % Store the L-test theta_bar as a warm-start for Phase 2
        registry.initial_ratios{1} = theta_bar_direct;
    else
        direct_estimation = false;

        % ---- Phase 1B: Binary search for t_star ----
        if params.verb >= 1
            fprintf('  Binary search for t* in [%.4f, %.4f]...\n', c0, t_max);
        end

        [t_star, ~, r_init_star] = binary_search_t(sampler_full, c, c0, t_max, params, points);
        r_init_star
        if params.verb >= 1
            fprintf('    t* = %.4f\n', t_star);
        end

        % Create body at t_star
        if t_star < t_max - 1e-8 * abs(t_max)
            registry = registry_store(registry, t_star, c, problem, ...
                                      polytope_opts, sampler_opts, sampler_full);
            % Store warm-start ratio for body P → t_star
            r_init_star
            if ~isempty(r_init_star)
                registry.initial_ratios{1} = r_init_star;
            end
            sampler_curr = registry.samplers{end};
            t_curr = t_star;
            if params.verb >= 1
                fprintf('    Created Polytope + Sampler for P ∩ H(t*).\n');
            end
        else
            % t_star ≈ t_max: use the full P sampler, skip intermediate bodies
            sampler_curr = sampler_full;
            t_curr = t_max;
            % Store the r_init for the full-P → P0 ratio (if available)
            if ~isempty(r_init_star)
                registry.initial_ratios{1} = r_init_star;
            end
            if params.verb >= 1
                fprintf('    t* ≈ t_max, using full P sampler.\n');
            end
        end
    end

    % ---- Phase 1C-D: Iteratively build intermediate bodies ----
    % Only run when we have an intermediate body and the L-test to P0
    % did not already pass on the initial batch.
    if ~direct_estimation && exist('sampler_curr', 'var')
        if params.verb >= 1
            fprintf('  Iterative construction (interval [c0=%.4f, t_curr=%.4f])...\n', ...
                    c0, t_curr);
        end

        for iter = 1:params.max_annealing_iter
            % Step D: check if we can reach P_0 directly
            [pts_check, sampler_curr] = advance_sampler(sampler_curr, n_needed);
            [~, L_pass, theta_bar, ~] = ul_test(pts_check, c0, c, params);

            if params.verb >= 2
                fprintf('    Iter %d: L-test to P0: theta_bar=%.4f, L_pass=%d\n', ...
                        iter, theta_bar, L_pass);
            end

            if L_pass
                if params.verb >= 1
                    fprintf('    L-test to target body passed. Schedule complete.\n');
                end
                registry.initial_ratios{end+1} = theta_bar;
                break;
            end

            % Binary search for next threshold in [c0, t_curr]
            [t_next, last_body, r_init_next] = binary_search_t(sampler_curr, c, c0, t_curr, params, pts_check);

            if params.verb >= 1
                fprintf('    Iter %d: t_next = %.4f\n', iter, t_next);
            end

            if last_body
                if params.verb >= 1
                    fprintf('    Last body reached; stopping phase construction.\n');
                end
                if ~isempty(r_init_next)
                    % registry.initial_ratios{k} is the ratio from body k to body k+1
                    r_init_next
                    registry.initial_ratios{end+1} = r_init_next;
                end
                break;
            end

            % Stop if t_next converged to an interval boundary
            %if t_next >= t_curr - 1e-8 * abs(t_curr) || t_next <= c0 + 1e-8
            %    if params.verb >= 1
            %        fprintf('    t_next converged to boundary; stopping iteration.\n');
            %    end
            %    break;
            %end

            % Create body at t_next
            registry = registry_store(registry, t_next, c, problem, ...
                                      polytope_opts, sampler_opts, sampler_curr);
            % Store warm-start ratio for this body pair
            if ~isempty(r_init_next)
                % registry.initial_ratios{k} is the ratio from body k to body k+1
                r_init_next
                registry.initial_ratios{end+1} = r_init_next;
            end
            sampler_curr = registry.samplers{end};
            t_curr = t_next;
        end

        if iter >= params.max_annealing_iter && ~L_pass
            warning('Annealing reached max iterations (%d) without L-test passing.', ...
                    params.max_annealing_iter);
        end
    elseif ~direct_estimation
        % t_star ≈ t_max and no sampler_curr: this is a degenerate case.
        % Fall back to direct estimation from full P.
        if params.verb >= 1
            fprintf('  No intermediate bodies; estimating directly from P.\n');
        end
    end

    n_bodies = length(registry.samplers);
    n_ratios = n_bodies;  % one ratio per sampler (outer body)

    if params.verb >= 1
        fprintf('  Annealing complete: %d samplers, %d ratios to estimate.\n', ...
                n_bodies, n_ratios);
        fprintf('  Thresholds (descending): ');
        for i = 1:length(registry.thresholds)
            fprintf('%.4f ', registry.thresholds{i});
        end
        fprintf('\n');
    end

    % =====================================================================
    % Phase 2: Estimate each volume ratio with sliding window
    % =====================================================================
    if params.verb >= 1
        fprintf('\n=== Phase 2: Estimating Volume Ratios ===\n');
    end

    eps_i = params.epsilon / sqrt(n_ratios + 1);
    log_ratios = zeros(n_ratios, 1);
    ratios = zeros(n_ratios, 1);
    steps_per_phase = zeros(n_ratios, 1);
    converged_flags = false(n_ratios, 1);

    for pair_idx = 1:n_ratios
        % Outer body sampler
        sampler = registry.samplers{pair_idx};

        % Inner threshold
        if pair_idx < n_ratios
            t_inner = registry.thresholds{pair_idx + 1};
        else
            t_inner = c0;  % target body P_0
        end

        if params.verb >= 1
            fprintf('  Ratio %d/%d: outer threshold=%s, inner threshold=%.4f\n', ...
                    pair_idx, n_ratios, ...
                    format_threshold(registry.thresholds{pair_idx}), t_inner);
        end

        % Retrieve warm-start ratio from Phase 1 (if any)
        %registry.initial_ratios
        if pair_idx <= length(registry.initial_ratios)
            r_init = registry.initial_ratios{pair_idx};
        else
            r_init = [];
        end

        % Estimate ratio via rejection sampling with sliding window
        r_init
        [r_hat, converged, n_steps] = estimate_ratio_sliding(...
            sampler, c, t_inner, eps_i, params, r_init);

        ratios(pair_idx) = r_hat;
        log_ratios(pair_idx) = log(max(r_hat, 1e-300));
        steps_per_phase(pair_idx) = n_steps;
        converged_flags(pair_idx) = converged;

        if params.verb >= 1
            fprintf('    r_hat = %.4f, converged = %d, steps = %d\n', ...
                    r_hat, converged, n_steps);
        end

        % Early exit if ratio is numerically zero
        if r_hat < 1e-300
            if params.verb >= 1
                fprintf('    Ratio numerically zero; stopping Phase 2.\n');
            end
            % Zero out remaining ratios
            log_ratios(pair_idx+1:end) = 0;
            ratios(pair_idx+1:end) = 0;
            converged_flags(pair_idx+1:end) = false;
            break;
        end
    end

    % =====================================================================
    % Phase 3: Assemble result
    % =====================================================================
    log_rho = sum(log_ratios);

    if log_rho < log(realmin) + 10
        warning('log(rho) = %.2f is near underflow; returning rho = 0.', log_rho);
        rho = 0.0;
    else
        rho = exp(log_rho);
    end

    result = build_result(registry, ratios, log_ratios, steps_per_phase, ...
                          n_ratios, log_rho, converged_flags);

    if params.verb >= 1
        fprintf('\n=== Complete ===\n');
        fprintf('  rho_hat   = %.16f\n', rho);
        fprintf('  log_rho   = %.16f\n', log_rho);
        fprintf('  n_phases  = %d\n', n_ratios);
        fprintf('  total HMC steps = %d\n', sum(steps_per_phase));
    end
end

% =========================================================================
% Subfunctions
% =========================================================================

function [polytope_opts, sampler_opts] = build_opts(params)
    %BUILD_OPTS Construct Polytope and Sampler options structs.
    %
    % Both Polytope() and Sampler() constructors accept the full options
    % struct and silently ignore fields they don't consume.  We build one
    % base struct from default_options() and fork it for each constructor.

    base_opts = default_options();

    % ---- Polytope options ----
    polytope_opts = base_opts;
    polytope_opts.presolve.runSimplify = true;
    polytope_opts.presolve.logFunc = @(tag, msg) 0;   % suppress constructor output

    % ---- Sampler options ----
    sampler_opts = base_opts;
    sampler_opts.maxTime = 86400 * 7;     % 1 week (prevent timeout)
    sampler_opts.maxStep = Inf;           % never stop on step count
    sampler_opts.N = Inf;                 % never stop on ESS
    sampler_opts.simdLen = params.simdLen;
    sampler_opts.initalStepSize = 0.2;
    sampler_opts.freezeMCMCAfterSamples = Inf;

    % Memory limit for chain storage (Module: MemoryStorage)
    sampler_opts.MemoryStorage.memoryLimit = 32 * 1024^3;  % 32 GB

    % Remove ProgressBar (cosmetic) but keep MixingTimeEstimator and MemoryStorage
    sampler_opts.module = setdiff(sampler_opts.module, {'ProgressBar'}, 'stable');

    % Suppress all sampler logging.  We control convergence ourselves via
    % the sliding window — the built-in ESS/step logging is not needed.
    % Note: this also suppresses diagnostic messages from DynamicStepSize
    % (e.g. step-shrink events).  Set params.verb >= 3 in the future to
    % route those through a function-handle logger.
    sampler_opts.logging = [];
end

function registry = registry_store(registry, t, c, problem, polytope_opts, sampler_opts, ...
                                   prev_sampler)
    %REGISTRY_STORE Create Polytope + Sampler for P ∩ H(t) and store in registry.
    %
    % If t == Inf, creates the full polytope P without adding a halfspace.
    % Polytope and Sampler are handle objects — mutations persist through the registry.
    %
    % Optional prev_sampler — if provided, its adapted stepSize is used as the
    % initial step size for the new sampler (warm-start across nested bodies).

    if isinf(t) && t > 0
        % Full polytope P (no halfspace constraint)
        problem_new = problem;
    else
        % Add halfspace c'*x <= t
        problem_new = problem;
        c_row = c(:)';  % ensure row vector for Aineq
        if ~isfield(problem_new, 'Aineq') || isempty(problem_new.Aineq)
            problem_new.Aineq = c_row;
            problem_new.bineq = t;
        else
            problem_new.Aineq = [problem_new.Aineq; c_row];
            problem_new.bineq = [problem_new.bineq(:); t];
        end
    end

    % Set a random seed for each new body (avoid identical chains across bodies)
    % Use uint32 to match the SimdTwister RNG expected by the Sampler constructor
    sampler_opts.seed = randi(intmax('uint32'), 'uint32');

    % Carry forward adapted step size from previous (outer) body.
    % Must be set before Sampler() reads opts.initalStepSize.
    if nargin >= 8 && ~isempty(prev_sampler) && isvalid(prev_sampler)
        sampler_opts.initalStepSize = prev_sampler.stepSize;
    end

    % PolytopeSamplerMatlab's top-level sample() initializes a shared timeout
    % clock before constructing Polytope and Sampler.  This routine constructs
    % those objects directly, so reproduce that required interface setup here.
    start_time = tic;
    polytope_opts.startTime = start_time;
    sampler_opts.startTime = start_time;

    polytope = Polytope(problem_new, polytope_opts);
    % sample() also initializes the SIMD random stream before Sampler().
    rng(sampler_opts.seed, 'simdTwister');
    sampler  = Sampler(polytope, sampler_opts);

    % ---- Burn-in: settle step size and approach stationarity, then freeze ----
    % Kook et al. (2022) report empirically O(d^0.52) steps per effective
    % sample for CRHMC on metabolic polytopes.  15 effective samples gives
    % the chain time to mix and DynamicStepSize time to stabilise.
    tic
    d = size(sampler.x, 2);
    burn_in_steps = ceil(15 * d^0.52);
    for k = 1:burn_in_steps
        sampler.step();
        if sampler.terminate == 3
            error('Sampler terminated with error during burn-in at step %d.', k);
        elseif sampler.terminate ~= 0
            warning('Sampler terminated with code %d during burn-in at step %d.', ...
                    sampler.terminate, k);
            break;
        end
    end
    time_to_burn_in = toc
    sampler.freezed = true;

    registry.thresholds{end+1} = t;
    registry.polytopes{end+1}  = polytope;
    registry.samplers{end+1}   = sampler;

    if registry.dim == 0
        registry.dim = size(sampler.x, 2);
    end
end

function [points, sampler] = advance_sampler(sampler, n_points)
    %ADVANCE_SAMPLER Advance the sampler and collect transformed points.
    %
    % Input:
    %   sampler  - Sampler handle object
    %   n_points - number of points to collect
    %
    % Output:
    %   points  - M x d_orig matrix of points in ORIGINAL space
    %             where M = ceil(n_points / simdLen) * simdLen
    %   sampler - the same handle object (mutated in place)

    simdLen = size(sampler.x, 1);
    if isempty(simdLen) || simdLen == 0
        simdLen = 1;
    end
    n_steps = ceil(n_points / simdLen);
    n_total = n_steps * simdLen;

    % Pre-allocate in original space
    dim_orig = size(sampler.problem.T, 1);
    points = zeros(n_total, dim_orig);

    idx = 0;  % last index written

    for step = 1:n_steps
        sampler.step();

        % Check for sampler errors
        if sampler.terminate == 3
            error('Sampler terminated with error at iteration %d.', sampler.i);
        elseif sampler.terminate ~= 0
            warning('Sampler terminated with code %d at iteration %d.', ...
                    sampler.terminate, sampler.i);
            break;
        end

        % Transform each chain's point to original space
        %   x_orig = T * x_internal + y
        for chain = 1:simdLen
            idx = (step-1)*simdLen + chain;
            x_int = sampler.x(chain, :)';
            points(idx, :) = (sampler.problem.T * x_int + sampler.problem.y)';
        end
    end

    % Trim if we collected fewer than allocated
    if idx < n_total
        points = points(1:idx, :);
    end
end

function [t_found, last_body, r_init] = binary_search_t(sampler, c, t_lo, t_hi, params, initial_points)
    %BINARY_SEARCH_T Find threshold t such that vol ratio ∈ [r, r+δ].
    %
    % Searches in [t_lo, t_hi] using U/L tests on points from the sampler.
    % The sampler is already sampling from the OUTER body.
    %
    % Optional: initial_points (M×d matrix) — if provided, used for the
    % first bisection iteration instead of collecting fresh points.
    %
    % Output:
    %   t_found    - threshold satisfying U-test and L-test
    %   last_body  - true if the L-test to t_lo passes (no more bodies needed)
    %   r_init     - theta_bar at t_found (warm-start for Phase 2)

    last_body = false;
    r_init = [];         % warm-start ratio for Phase 2
    lo = t_lo;
    hi = t_hi;
    t_found = (lo + hi) / 2;  % fallback
    have_initial = (nargin >= 6 && ~isempty(initial_points));
    % Use pre-collected points for first iteration, otherwise collect fresh
    if have_initial
        pts = initial_points;
    else
        n_needed = params.nu * params.N_utest;
        [pts, ~] = advance_sampler(sampler, n_needed);
    end
    %[~, L_pass, theta_bar_lo, ~] = ul_test(pts, t_lo, c, params);
    %if L_pass
    %    t_found = t_lo;
    %    last_body = true;
    %    r_init = theta_bar_lo;
    %    return;
    %end
    theta_bar = 0;
    for i = 1:params.max_bisect_iter
        if i==params.max_bisect_iter
            fprintf('   hello world max_bisect_iter reached\n');
        end
        if hi - lo < 1e-8 * max(abs(hi), 1)
            t_found = (lo + hi) / 2;
            fprintf('   hello world t_found = %.4f\n', t_found);
            r_init = theta_bar;
            break;
        end

        t_mid = (lo + hi) / 2;

        [U_pass, L_pass, theta_bar, s_theta] = ul_test(pts, t_mid, c, params);

        if params.verb >= 2
            fprintf('      bisect %d: t=%.4f, theta_bar=%.4f, s=%.4f, U=%d, L=%d\n', ...
                    i, t_mid, theta_bar, s_theta, U_pass, L_pass);
        end

        if U_pass && L_pass
            t_found = t_mid;
            r_init = theta_bar;
            return;
        elseif ~U_pass
            % Ratio >= r+delta → inner body too large → tighten (decrease t)
            hi = t_mid;
        else  % ~L_pass
            % Ratio <= r → inner body too small → loosen (increase t)
            lo = t_mid;
        end
    end

    % If we exit without exact convergence, use midpoint
    if ~exist('U_pass', 'var') || ~(U_pass && L_pass)
        t_found = (lo + hi) / 2;
    end
end

function [r_hat, converged, n_steps] = estimate_ratio_sliding(sampler, c, t_inner, ...
                                                              eps_i, params, r_init)
    %ESTIMATE_RATIO_SLIDING Estimate volume ratio via rejection sampling with sliding window.
    %
    % Samples from the outer body (sampler) and tests membership in the
    % inner body (c'*x <= t_inner). Convergence via sliding window.
    %
    % Optional r_init (scalar in (0,1)) warm-starts the counters and
    % pre-fills the sliding window from a Phase 1 U/L-test ratio estimate.

    % Use this sampler's internal dimension for window sizing
    dim = size(sampler.x, 2);
    window = new_window(params, dim, eps_i);
    simdLen = size(sampler.x, 1);
    if isempty(simdLen) || simdLen == 0
        simdLen = 1;
    end
    c_vec = c(:);

    % ---- Warm-start from Phase 1 ratio estimate ----
    if nargin >= 7 && ~isempty(r_init) && r_init > 0 && r_init < 1
        % Use the same sample size as Phase 1 for consistent weighting
        n_warm = params.nu * params.N_utest;
        count_in = round(r_init * n_warm);
        total = n_warm;
        % Pre-fill the sliding window with the initial estimate
        %window.last_W(:) = r_init;
        %window.min_val = r_init;
        %window.max_val = r_init;
        %window.min_index = 1;
        %window.max_index = 1;
        %window.index = 1;
        if params.verb >= 2
            fprintf('      (warm-start: r_init=%.4f, n_warm=%d)\n', r_init, n_warm);
        end
    else
        count_in = 0;
        total = 0;
    end
    n_steps = 0;
    converged = false;

    while ~converged
        sampler.step();
        n_steps = n_steps + simdLen;

        % Check for errors
        if sampler.terminate == 3
            error('Sampler terminated with error at iteration %d.', sampler.i);
        elseif sampler.terminate ~= 0
            warning('Sampler terminated with code %d at iteration %d.', ...
                    sampler.terminate, sampler.i);
            break;
        end

        for chain = 1:simdLen
            % Transform from internal to original space
            x_int = sampler.x(chain, :)';
            x_orig = sampler.problem.T * x_int + sampler.problem.y;

            % Membership test
            in_smaller = (c_vec' * x_orig <= t_inner);

            count_in = count_in + in_smaller;
            total = total + 1;

            % Update running average and sliding window
            val = count_in / total;
            window = update_window(window, val);

            if window.converged
                converged = true;
                break;
            end
        end
    end

    r_hat = count_in / total;
end

function result = build_result(registry, ratios, log_ratios, steps_per_phase, ...
                               n_ratios, log_rho, converged)
    %BUILD_RESULT Assemble the output result struct.
    result.thresholds = registry.thresholds;
    result.ratios = ratios;
    result.log_ratios = log_ratios;
    result.steps_per_phase = steps_per_phase;
    result.total_steps = sum(steps_per_phase);
    result.n_phases = n_ratios;
    result.log_rho = log_rho;
    result.converged = converged;
    result.samplers = registry.samplers;
end

function s = format_threshold(t)
    %FORMAT_THRESHOLD String representation of a threshold value.
    if isinf(t) && t > 0
        s = '+Inf';
    else
        s = sprintf('%.4f', t);
    end
end
