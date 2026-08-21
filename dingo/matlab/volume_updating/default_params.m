function params = default_params()
    %DEFAULT_PARAMS Return the default parameter struct for relative_volume.
    %
    % Usage:
    %   params = default_params();
    %   params = Setfield(default_params(), user_overrides);

    % ---- Annealing schedule parameters ----
    params.r     = 0.1;   % target lower bound for each volume ratio
    params.delta = 0.05;   % half-width of target ratio interval [r, r+delta]
    params.alpha = 0.20;   % significance level for U/L t-tests
    params.nu    = 10;     % number of sublists for U/L tests
    params.N_utest = 320;  % points per sublist (total = nu * N_utest)
    params.window_size = 5000;  % size of the sliding window (We now set it in the new_window()!!)

    % ---- Convergence parameters ----
    params.epsilon = 0.10;             % target relative error for final estimate
    params.max_annealing_iter = 200;    % max Phase 1 construction iterations
    params.max_bisect_iter = 100;       % max binary search iterations per body

    % ---- Sampler parameters ----
    params.simdLen = 8;                % number of parallel chains per step()

    % ---- t_max estimation ----
    params.t_max_safety = 1.01;        % inflate sampled t_max by this factor

    % ---- Verbosity ----
    params.verb = 1;   % 0 = quiet, 1 = normal, 2 = detailed
end
