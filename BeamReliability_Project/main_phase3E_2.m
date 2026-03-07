clear;
clc;
addpath(genpath(pwd));

fprintf('==============================\n')
fprintf('PHASE 3E_2 — MASTER RELIABILITY DIAGRAM\n')
fprintf('==============================\n')

% Base parameters
params.mu_E = 210e9;
params.mu_F = 1000;

params.sigma_E = 0.05 * params.mu_E;
params.sigma_F = 0.05 * params.mu_F;

params.L = 2;
params.I = 8e-6;
params.delta_allow = 0.005;

params.verbose = false;

% ---------------------------------
% Master Reliability Phase Diagram
% ---------------------------------
fprintf('\nGenerating reliability phase diagram...\n')
master_reliability_phase_diagram(params);

% ---------------------------------
% Sampling Efficiency Comparison
% ---------------------------------
fprintf('\nRunning sampling efficiency comparison...\n')
compare_sampling_efficiency(params);

% ---------------------------------
% Reliability atlas
% ---------------------------------
fprintf('\nGenerating reliability atlas...\n')
reliability_atlas(params);

% ---------------------------------
% Reliability geometry diagram
% ---------------------------------
fprintf('\nGenerating reliability geometry diagram...\n')
reliability_geometry_diagram(params);

fprintf('\nPhase 3E_2 completed.\n')