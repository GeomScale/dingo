function simdLen = init()
    %INIT Add PolytopeSamplerMatlab and src/ to the MATLAB path and compile solvers.
    %
    %   simdLen = init()
    %
    % Returns the maximum usable simdLen (number of parallel SIMD chains).
    % Tries simdLen=8 first; falls back to simdLen=4 if the MEX solver for 8
    % cannot be compiled (e.g. no C++ compiler is available).  simdLen=4 has
    % pre-compiled binaries in PolytopeSamplerMatlab/bin/ for Linux, macOS,
    % and Windows on x86-64.
    %
    % This function is idempotent — subsequent calls return the cached value
    % without recompiling.  Called automatically by relative_volume() and
    % run_relative_volume().

    persistent cached_simdLen;
    if ~isempty(cached_simdLen)
        simdLen = cached_simdLen;
        return;
    end

    src_dir = fileparts(mfilename('fullpath'));
    root    = fileparts(src_dir);

    % Add src/ so default_params.m, sliding_window.m, ul_test.m are findable
    addpath(src_dir);

    % The original repository keeps PolytopeSamplerMatlab beside src/.  The
    % Dingo overlay can be installed in a different layout, so first honor an
    % explicit environment variable and then check deterministic locations.
    sampler_root = getenv('DINGO_POLYTOPE_SAMPLER_MATLAB');
    candidates = {
        sampler_root, ...
        fullfile(root, 'PolytopeSamplerMatlab'), ...
        fullfile(fileparts(fileparts(fileparts(src_dir))), 'PolytopeSamplerMatlab'), ...
        fullfile(fileparts(fileparts(fileparts(fileparts(src_dir)))), 'PolytopeSamplerMatlab') ...
    };
    sampler_root = '';
    for candidate_idx = 1:numel(candidates)
        candidate = candidates{candidate_idx};
        if ~isempty(candidate) && exist(fullfile(candidate, 'code'), 'dir')
            sampler_root = candidate;
            break;
        end
    end
    if isempty(sampler_root)
        error(['PolytopeSamplerMatlab was not found. Set ', ...
               'DINGO_POLYTOPE_SAMPLER_MATLAB to its repository root.']);
    end
    sampler_code = fullfile(sampler_root, 'code');
    sampler_bin  = fullfile(sampler_root, 'bin');

    addpath(genpath(sampler_code));
    addpath(sampler_bin);

    % Scalar solver — pre-compiled binaries exist for all platforms
    compile_solver(0);

    % SIMD solver.  Try simdLen=8 first (best throughput); fall back to
    % simdLen=4 if the MEX binary cannot be found or compiled.  simdLen=4
    % has pre-compiled binaries shipped with PolytopeSamplerMatlab.
    for trial = [8, 4]
        try
            compile_solver(trial);
            cached_simdLen = trial;
            if trial == 4
                warning('init:fallback', ...
                        ['Cannot compile MEX solver for simdLen=8 ' ...
                         '(no C++ compiler?).  Falling back to simdLen=4 ' ...
                         'which uses pre-compiled binaries.']);
            end
            simdLen = cached_simdLen;
            return;
        catch ME
            if trial == 4
                rethrow(ME);
            end
            % simdLen=8 failed — try simdLen=4
        end
    end
end
