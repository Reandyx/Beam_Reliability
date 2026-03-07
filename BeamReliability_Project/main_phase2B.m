%% ==============================================================
%  PHASE 2B — FORM (2D Reliability, HL-RF)
%  Random Variables: E ~ N(mu_E, sigma_E)
%                    F ~ N(mu_F, sigma_F)
%  Method: Hasofer-Lind-Rackwitz-Fiessler (HL-RF)
% ==============================================================

clear; clc;
addpath(genpath(pwd));

%% --------------------------------------------------------------
% 1. Deterministic Beam Parameters
% --------------------------------------------------------------

L = 2.0;             % [m]
b = 0.05;            % [m]
h = 0.10;            % [m]
delta_allow = 0.0032;

I = section_inertia_rect(b,h);

%% --------------------------------------------------------------
% 2. Stochastic Parameters
% --------------------------------------------------------------

mu_E = 210e9;
sigma_E = 0.05 * mu_E;

mu_F = 1000;
sigma_F = 0.10 * mu_F;

% Pack into parameter struct
params.L = L;
params.I = I;
params.delta_allow = delta_allow;

params.mu_E = mu_E;
params.sigma_E = sigma_E;

params.mu_F = mu_F;
params.sigma_F = sigma_F;
params.F = mu_F;   % Deterministic F for 1D solver compatibility

%% --------------------------------------------------------------
% 3. 2D Monte Carlo Simulation
% --------------------------------------------------------------

N = 1e6;

[E_samples, F_samples] = generate_EF_samples( ...
    mu_E, sigma_E, mu_F, sigma_F, N);

[Pf_MC, CI_MC] = failure_probability_mc_2d( ...
    E_samples, F_samples, params);

%% --------------------------------------------------------------
% 4. FORM — HL-RF
% --------------------------------------------------------------

[beta_FORM, Z_star, Pf_FORM] = form_2d_hlrf(params);

alpha = importance_factors(Z_star, beta_FORM);

%% --------------------------------------------------------------
% 5. Reduction Consistency Test (σF → small)
% --------------------------------------------------------------

params_test = params;
params_test.sigma_F = 0.01 * params.mu_F;  % 1% CoV

[beta_2D_small, ~, ~] = form_2d_hlrf(params_test);
[beta_1D, ~, ~, ~, ~] = form_1d_solver(params);

reduction_error = abs(beta_2D_small - beta_1D);
rel_reduction_error = reduction_error / beta_1D;

%% --------------------------------------------------------------
% 6. Structured Output
% --------------------------------------------------------------

fprintf('\n=========================================================\n');
fprintf('PHASE 2B — FORM (2D HL-RF Implementation)\n');
fprintf('=========================================================\n');

fprintf('\n--- Monte Carlo (2D) ---\n');
fprintf('Pf_MC         = %.6e\n', Pf_MC);
fprintf('CI (95%%)      = ± %.6e\n', CI_MC);

fprintf('\n--- FORM (HL-RF) ---\n');
fprintf('Beta           = %.6f\n', beta_FORM);
fprintf('Pf_FORM        = %.6e\n', Pf_FORM);
fprintf('Design point Z* = [%.6f , %.6f]\n', Z_star(1), Z_star(2));

fprintf('\n--- Importance Factors ---\n');
fprintf('alpha_E        = %.6f\n', alpha(1));
fprintf('alpha_F        = %.6f\n', alpha(2));

if abs(alpha(1)) > abs(alpha(2))
    fprintf('Dominant source: Material variability (E)\n');
else
    fprintf('Dominant source: Load variability (F)\n');
end

fprintf('\n--- Validation ---\n');

fprintf('Reduction error (abs)      = %.6e\n', reduction_error);
fprintf('Reduction error (relative) = %.6e\n', rel_reduction_error);

if rel_reduction_error < 0.03
    fprintf('✓ 2D formulation reduces correctly to 1D (σF → small)\n');
else
    fprintf('⚠ Reduction discrepancy detected\n');
end

abs_Pf_error = abs(Pf_MC - Pf_FORM);
fprintf('MC vs FORM Pf difference   = %.6e\n', abs_Pf_error);

if abs_Pf_error < 3*CI_MC
    fprintf('✓ MC and FORM agree within statistical tolerance\n');
else
    fprintf('⚠ Deviation noticeable (possible nonlinearity effects)\n');
end

fprintf('=========================================================\n');