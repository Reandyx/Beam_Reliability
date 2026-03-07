%% ============================================================
%  PHASE 2A — FORM (1D Reliability)
%  Monte Carlo vs FORM Comparison
%  Single Random Variable: E ~ N(mu_E, sigma_E)
% =============================================================

clear; clc;

addpath(genpath(pwd));

%% -------------------------------------------------------------
% 1️⃣ Deterministic Beam Parameters
%% -------------------------------------------------------------

F = 1000;            % Load [N]
L = 2.0;             % Length [m]
b = 0.05;            % Width [m]
h = 0.10;            % Height [m]
delta_allow = 0.0032;  % adjusted to position system near limit state for reliability analysis

I = section_inertia_rect(b,h);

%% -------------------------------------------------------------
% 2️⃣ Stochastic Parameters
%% -------------------------------------------------------------

mu_E    = 210e9;      % Mean Young modulus [Pa]
sigma_E = 0.05*mu_E;  % 5% CoV

% ---- Parameter struct (Interface hardening) ----
params.F = F;
params.L = L;
params.I = I;
params.delta_allow = delta_allow;
params.mu_E = mu_E;
params.sigma_E = sigma_E;

%% -------------------------------------------------------------
% 3️⃣ Monte Carlo Simulation
%% -------------------------------------------------------------

N = 1e6;

% Generate Young modulus samples
E_samples = generate_E_samples(mu_E, sigma_E, N);

% Compute deflections
delta_samples = monte_carlo_deflection(F, L, E_samples, I);

% Compute failure probability
[Pf_MC, CI_MC] = failure_probability_mc(delta_samples, delta_allow);

%% -------------------------------------------------------------
% 4️⃣ FORM Computation
%% -------------------------------------------------------------

[beta_FORM, beta_signed, Pf_FORM, Z_star, E_star] = ...
    form_1d_solver(params);

beta_analytical = reliability_index(params);

%% -------------------------------------------------------------
% 5️⃣ Validation
%% -------------------------------------------------------------

abs_beta_error = abs(beta_FORM - beta_analytical);
rel_beta_error = abs_beta_error / beta_analytical;

abs_Pf_error = abs(Pf_MC - Pf_FORM);
rel_Pf_error = abs_Pf_error / Pf_MC;

%% -------------------------------------------------------------
% 6️⃣ Structured Output (Publication Order)
%% -------------------------------------------------------------

fprintf('\n============================================================\n');
fprintf('PHASE 2A — FORM (1D Reliability Implementation)\n');
fprintf('============================================================\n');

%% --- Deterministic Parameters ---
fprintf('\n--- DETERMINISTIC PARAMETERS ---\n');
fprintf('F              = %.3f N\n', params.F);
fprintf('L              = %.3f m\n', params.L);
fprintf('I              = %.6e m^4\n', params.I);
fprintf('delta_allow    = %.6e m\n', params.delta_allow);

%% --- Stochastic Parameters ---
fprintf('\n--- STOCHASTIC PARAMETERS ---\n');
fprintf('mu_E           = %.6e Pa\n', params.mu_E);
fprintf('sigma_E        = %.6e Pa\n', params.sigma_E);
fprintf('CoV            = %.3f\n', params.sigma_E/params.mu_E);

%% --- Monte Carlo Results ---
fprintf('\n--- MONTE CARLO RESULTS ---\n');
fprintf('Pf_MC          = %.6e\n', Pf_MC);
fprintf('Abs Pf error   = %.6e\n', abs_Pf_error);
fprintf('Rel Pf error   = %.6e\n', rel_Pf_error);

%% --- FORM Results ---
fprintf('\n--- FORM RESULTS ---\n');
fprintf('Design point E*      = %.6e Pa\n', E_star);
fprintf('Design point Z*      = %.6f\n', Z_star);
fprintf('Beta (signed)        = %.6f\n', beta_signed);
fprintf('Beta (magnitude)     = %.6f\n', beta_FORM);
fprintf('Pf_FORM              = %.6e\n', Pf_FORM);

%% --- Validation ---
fprintf('\n--- VALIDATION ---\n');
fprintf('Beta abs diff        = %.6e\n', abs_beta_error);
fprintf('Beta rel diff        = %.6e\n', rel_beta_error);

if rel_beta_error < 1e-6
    fprintf('✔ Analytical validation PASSED\n');
else
    fprintf('✖ Analytical validation FAILED\n');
end

if abs_Pf_error < 5e-3*Pf_MC
    fprintf('✔ Monte Carlo vs FORM agreement within statistical tolerance\n');
else
    fprintf('⚠ Monte Carlo deviation noticeable (possible nonlinearity)\n');
end

%% --- Interpretation ---
fprintf('\n--- INTERPRETATION ---\n');

if beta_signed < 0
    fprintf('Failure governed by decrease in E (left tail).\n');
else
    fprintf('Failure governed by increase in E (right tail).\n');
end

if rel_Pf_error > 0.01
    fprintf('FORM linearization limits becoming noticeable.\n');
else
    fprintf('FORM approximation consistent with Monte Carlo.\n');
end

fprintf('\n============================================================\n');