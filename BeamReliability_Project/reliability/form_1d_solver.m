function [beta, beta_signed, Pf_FORM, Z_star, E_star] = form_1d_solver(params)

% Extract stochastic parameters
mu_E    = params.mu_E;
sigma_E = params.sigma_E;

% Transform limit state to standard normal space
g_Z = @(Z) limit_state_E(mu_E + sigma_E .* Z, params);

% ---- Root existence check ----
g_left  = g_Z(-5);
g_right = g_Z(5);

if sign(g_left) == sign(g_right)
    warning('No failure surface detected in search interval [-5,5].');
end

% ---- Robust root finding ----
Z_star = fzero(g_Z, [-5 5]);

% Reliability index
beta_signed = Z_star;
beta = abs(Z_star);

% Failure probability
Pf_FORM = normcdf(-beta);

% Back-transform to physical space
E_star = mu_E + sigma_E * Z_star;

% ---- Physical validation ----
g_check = limit_state_E(E_star, params);
if abs(g_check) > 1e-10
    warning('Design point validation failed: g(E*) not near zero.');
end

end