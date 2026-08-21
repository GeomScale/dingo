function [U_pass, L_pass, theta_bar, s_theta] = ul_test(points, t_inner, c, params)
    %UL_TEST U-test and L-test for volume ratio interval verification.
    %
    % Input:
    %   points  - M x d matrix of points in ORIGINAL space from the OUTER body
    %   t_inner - threshold defining the INNER body (c'*x <= t_inner)
    %   c       - cost vector (d x 1)
    %   params  - struct with fields: .r, .delta, .alpha, .nu, .N_utest
    %
    % Output:
    %   U_pass    - true if U-test rejects H0: theta >= r+delta
    %   L_pass    - true if L-test rejects H0: theta <= r
    %   theta_bar - sample mean of per-sublist empirical ratios
    %   s_theta   - sample standard deviation of per-sublist ratios
    %
    % Reference: Chalkis et al. (2023), Section 4.1, Eqs. 9-10.

    nu = params.nu;
    N  = params.N_utest;
    r  = params.r;
    delta = params.delta;
    alpha = params.alpha;

    n_available = size(points, 1);
    n_needed = nu * N;

    if n_available < n_needed
        warning('ul_test: only %d points available, need %d. Returning not-passed.', ...
                n_available, n_needed);
        U_pass = false;
        L_pass = false;
        theta_bar = NaN;
        s_theta = NaN;
        return;
    end

    % Compute c'*x for the first nu*N points
    c_vals = points(1:n_needed, :) * c(:);

    % Partition into nu sublists and compute per-sublist ratios
    theta_hat = zeros(nu, 1);
    for j = 1:nu
        idx_start = (j-1)*N + 1;
        idx_end   = j*N;
        n_in = sum(c_vals(idx_start:idx_end) <= t_inner);
        theta_hat(j) = n_in / N;
    end

    % Sample statistics
    theta_bar = mean(theta_hat);
    s_theta   = std(theta_hat);

    % Critical value: upper alpha-quantile of t-distribution with nu-1 df
    t_crit = tinv(1 - alpha, nu - 1);

    se = s_theta / sqrt(nu);   % standard error

    % U-test: reject H0: theta >= r+delta  if  theta_bar <= r+delta - t_crit*se
    U_pass = (theta_bar <= r + delta - t_crit * se);

    % L-test: reject H0: theta <= r         if  theta_bar >= r + t_crit*se
    L_pass = (theta_bar >= r + t_crit * se);
end
