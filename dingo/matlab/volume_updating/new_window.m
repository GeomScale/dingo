function window = new_window(params, dim, eps_i)
    %NEW_WINDOW Create a sliding window convergence tracker.
    %
    % Input:
    %   dim   - effective dimension of the polytope
    %   eps_i - per-ratio error budget
    %
    % Output:
    %   window - struct with fields:
    %     .W          window size = ceil(4*dim^2 + 500)
    %     .min_val    current minimum in the window
    %     .max_val    current maximum in the window
    %     .min_index  index of min in circular buffer
    %     .max_index  index of max in circular buffer
    %     .last_W     circular buffer of length W
    %     .index      write pointer (1-based)
    %     .converged  logical flag
    %     .eps_i      per-ratio error budget
    %
    % Reference: Cousins & Vempala (2016), Section 3.3.
    % Adapted from Volume-and-Sampling/Volume.m lines 246-266.

    %window.W = ceil(4 * dim^2 + 500);
    window.W = ceil(params.simdLen * dim^0.52)+400;
    %window.W = params.window_size;
    window.min_val = -realmax;
    window.max_val = realmax;
    window.min_index = window.W;
    window.max_index = window.W;
    window.last_W = zeros(window.W, 1);
    window.index = 1;
    window.converged = false;
    window.eps_i = eps_i;
    window.filled = false;
end

function window = update_window(window, val)
    %UPDATE_WINDOW Push a new cumulative running average into the window.
    %
    % Input:
    %   window - window struct from new_window
    %   val    - current cumulative running average (scalar)
    %
    % Output:
    %   window - updated struct with .converged set to true if converged
    %
    % Maintains O(1) amortized min/max tracking via lazy eviction:
    % when the outgoing value was NOT the min or max, updates are O(1).
    % When the outgoing value WAS the min or max, a full scan recomputes.
    %
    % Adapted from Volume-and-Sampling/Volume.m lines 268-301.

    window.last_W(window.index) = val;

    % Update minimum
    if val <= window.min_val
        window.min_val = val;
        window.min_index = window.index;
    elseif window.min_index == window.index
        % The outgoing element was the min -- recompute
        window.filled = true;
        [window.min_val, window.min_index] = min(window.last_W);
    end

    % Update maximum
    if val >= window.max_val
        window.max_val = val;
        window.max_index = window.index;
    elseif window.max_index == window.index
        % The outgoing element was the max -- recompute
        window.filled = true;
        [window.max_val, window.max_index] = max(window.last_W);
    end

    % Check convergence: relative spread within eps_i/2
    if (window.max_val - window.min_val) / window.max_val <= window.eps_i / 2 && window.filled
    %if (max(window.last_W) - min(window.last_W)) / max(window.last_W) <= window.eps_i / 2
        window.converged = true;
    end
    % Advance circular buffer pointer
    window.index = mod(window.index, window.W) + 1;
end
