function [Pf, CI] = failure_probability_mc_2d(E_samples, F_samples, params)

% FAILURE_PROBABILITY_MC_2D
% Monte Carlo estimation of failure probability for the 2D reliability problem.
% Random variables:
%   E - Young's modulus
%   F - Load
% Limit state:
%   g(E,F) = delta_allow - delta(E,F)
% Failure condition:
%   g < 0
% Vectorized implementation used for computational efficiency.
% Equivalent to loop-based Monte Carlo evaluation.

N = length(E_samples);

% Compute beam deflection for all samples (vectorized)
delta = beam_deflection(F_samples, params.L, E_samples, params.I);

% Evaluate limit state function
g = params.delta_allow - delta;

% Failure indicator
failures = (g < 0);

% Failure probability estimate
Pf = mean(failures);

% Avoid numerical issues when Pf = 0
Pf = max(Pf, eps);

% 95% confidence interval (normal approximation)
CI = 1.96 * sqrt(Pf * (1 - Pf) / N);

end