clear
clc

fprintf('==============================\n')
fprintf('PHASE 3F — GLOBAL SENSITIVITY\n')
fprintf('==============================\n')

%% PARAMETERS

params.mu_E = 210e9;
params.sigma_E = 21e9;

params.mu_F = 1000;
params.sigma_F = 200;

params.L = 2;
params.I = 8e-6;

%% Sobol sampling

N = 1e5;

fprintf('\nRunning Sobol sensitivity analysis...\n')

[S_E,S_F,S_E_total,S_F_total] = sobol_indices(params,N);

%% Results

fprintf('\nSobol First Order Indices\n')
fprintf('S_E = %.3f\n',S_E)
fprintf('S_F = %.3f\n',S_F)

fprintf('\nSobol Total Order Indices\n')
fprintf('S_E_total = %.3f\n',S_E_total)
fprintf('S_F_total = %.3f\n',S_F_total)

fprintf('\nCheck: S_E + S_F = %.3f\n',S_E + S_F)

%% Visualization

plot_sobol_results(S_E,S_F,S_E_total,S_F_total)

fprintf('\nPhase 3F completed successfully.\n')