function [beta_analytical, Pf_analytical] = analytical_reliability_estimate(params)

% Analytical reliability estimate using ratio R = F/E

mu_E = params.mu_E;
sigma_E = params.sigma_E;

mu_F = params.mu_F;
sigma_F = params.sigma_F;

L = params.L;
I = params.I;
delta_allow = params.delta_allow;

% ---- Critical ratio ----
K = (3 * I * delta_allow) / (L^3);

% ---- Mean of ratio ----
mu_R = mu_F / mu_E;

% ---- Variance of ratio (first-order approximation) ----
sigma_R2 = (sigma_F / mu_E)^2 + ((mu_F * sigma_E) / (mu_E^2))^2;
sigma_R = sqrt(sigma_R2);

% ---- Reliability index ----
beta_analytical = (K - mu_R) / sigma_R;

% ---- Failure probability ----
Pf_analytical = normcdf(-beta_analytical);

end