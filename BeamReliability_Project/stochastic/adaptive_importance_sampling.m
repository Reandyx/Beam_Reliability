function [Pf_est, variance_est, mu_history, ESS] = adaptive_importance_sampling(params, Z_star, Sigma, N, max_iter)

% ADAPTIVE_IMPORTANCE_SAMPLING
% Iterative adaptive importance sampling for failure probability
%
% Inputs:
% params   : parameter struct
% Z_star   : FORM design point
% Sigma    : covariance matrix
% N        : samples per iteration
% max_iter : AIS iterations
%
% Outputs:
% Pf_est        : failure probability estimate
% variance_est  : estimator variance
% mu_history    : trajectory of sampling centers
% ESS           : effective sample size

dim = length(Z_star);

% Initial sampling mean
mu = Z_star(:);

% Store trajectory
mu_history = mu';

% Relaxation parameter (stability)
alpha = 1;

% Storage for estimator accumulation
Pf_sum = 0;
total_samples = 0;

% Storage for variance estimation
all_weighted_samples = [];

for k = 1:max_iter

    % ---------------------------
    % Generate samples
    % ---------------------------
    U = mvnrnd(mu', Sigma, N);

    % ---------------------------
    % Transform to physical variables
    % ---------------------------
    Z_E = U(:,1);
    Z_F = U(:,2);

    E = params.mu_E + params.sigma_E .* Z_E;
    F = params.mu_F + params.sigma_F .* Z_F;

    % ---------------------------
    % Evaluate limit state
    % ---------------------------
    delta = beam_deflection(F, params.L, E, params.I);
    g = params.delta_allow - delta;

    % ---------------------------
    % Failure indicator
    % ---------------------------
    I = (g < 0);
    num_fail = sum(I);
    fprintf('AIS iteration %d failure samples = %d\n', k, num_fail);

    % ---------------------------
    % Importance weights
    % ---------------------------
    phi_original = mvnpdf(U,[0 0],eye(dim));
    phi_shifted  = mvnpdf(U,mu',Sigma);

    w = phi_original ./ phi_shifted;

    weighted_samples = I .* w;

    % ---------------------------
    % Accumulate estimator
    % ---------------------------
    Pf_sum = Pf_sum + sum(weighted_samples);
    total_samples = total_samples + N;

    all_weighted_samples = [all_weighted_samples; weighted_samples];

    % ---------------------------
    % Update sampling center
    % ---------------------------
    if sum(weighted_samples) > 0

        mu_new = sum(weighted_samples .* U) ./ sum(weighted_samples);

        % Relaxed update for stability
        mu = (1 - alpha) * mu + alpha * mu_new';

    end

    mu_history = [mu_history; mu'];

end

% ---------------------------
% Final probability estimate
% ---------------------------
Pf_est = Pf_sum / total_samples;

% ---------------------------
% Variance estimate
% ---------------------------
variance_est = var(all_weighted_samples) / total_samples;

% ---------------------------
% Effective sample size
% ---------------------------
ESS = compute_ess(all_weighted_samples);

end