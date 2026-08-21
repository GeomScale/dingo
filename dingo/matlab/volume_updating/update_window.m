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
        [window.min_val, window.min_index] = min(window.last_W);
    end

    % Update maximum
    if val >= window.max_val
        window.max_val = val;
        window.max_index = window.index;
    elseif window.max_index == window.index
        % The outgoing element was the max -- recompute
        [window.max_val, window.max_index] = max(window.last_W);
    end

    % Check convergence: relative spread within eps_i/2
    if (window.max_val - window.min_val) / window.max_val <= window.eps_i / 2
        window.converged = true;
    end

    % Advance circular buffer pointer
    window.index = mod(window.index, window.W) + 1;
end
