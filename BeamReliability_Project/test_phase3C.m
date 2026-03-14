clear; clc;
addpath(genpath(pwd));

params.verbose = false; 

params.mu_E = 210e9;
params.sigma_E = 0.1 * params.mu_E;

params.mu_F = 3000;
params.sigma_F = 0.1 * params.mu_F;
params.F = params.mu_F;

params.L = 2;
params.I = 8.333e-6;

params.delta_allow = 0.005;

%% Analytical reliability
[beta_a,Pf_a] = analytical_reliability_estimate(params)

%% FORM
[beta_FORM,~,Pf_FORM] = form_1d_solver(params)

%% SORM
Pf_SORM = sorm_correction(params)

%% Method comparison plot
method_comparison_plot(params)

%% Scaling law verification
scaling_law_verification(params)

%% Reliability regime diagram
reliability_regime_diagram(params)