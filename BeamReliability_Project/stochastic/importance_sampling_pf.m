function [Pf_IS, variance_IS, ESS, U_samples, g_vals, weights] = ...
    importance_sampling_pf(params, Z_star, Sigma, N)

% IMPORTANCE_SAMPLING_PF
% Baseline importance sampling estimator for failure probability
%
% Inputs
% params  : parameter struct
% Z_star  : FORM design point
% Sigma   : sampling covariance matrix
% N       : number of samples
%
% Outputs
% Pf_IS        : failure probability estimate
% variance_IS  : estimator variance
% ESS          : effective sample size
% U_samples    : samples in standard space
% g_vals       : limit state values
% weights      : importance weights

dim = length(Z_star);

% Generate samples
U_samples = mvnrnd(Z_star', Sigma, N);

% Transform to physical variables
Z_E = U_samples(:,1);
Z_F = U_samples(:,2);

E = params.mu_E + params.sigma_E .* Z_E;
F = params.mu_F + params.sigma_F .* Z_F;

% Evaluate limit state
delta = beam_deflection(F, params.L, E, params.I);
g_vals = params.delta_allow - delta;

% Failure indicator
I = (g_vals < 0);
num_fail = sum(I);
fprintf('IS failure samples = %d\n', num_fail);

% Importance weights
phi_original = mvnpdf(U_samples,[0 0],eye(dim));
phi_shifted  = mvnpdf(U_samples,Z_star',Sigma);

weights = phi_original ./ phi_shifted;

% Probability estimate
Pf_IS = mean(I .* weights);

% Variance estimate
variance_IS = var(I .* weights) / N;

% Effective sample size
ESS = compute_ess(weights .* I);

end