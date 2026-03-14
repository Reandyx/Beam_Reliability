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

% --- Phase 3D analysis ---

[U_history, design_point, curvature] = hlrf_trajectory_analysis(params);

fprintf('Design Point:\n');
disp(design_point)

fprintf('Curvature measure: %.6e\n',curvature);

% --- Plot geometry ---

plot_hlrf_geometry(U_history,design_point,params);