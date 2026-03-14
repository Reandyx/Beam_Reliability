%% ============================================================
% PHASE 2C — RELIABILITY REGIME MAPPING
% Study reliability behavior under uncertainty regimes
%
% Random Variables:
%   E ~ N(mu_E, sigma_E)
%   F ~ N(mu_F, sigma_F)
%
% Grid:
%   CoV_E = 0.02 : 0.02 : 0.15
%   CoV_F = 0.02 : 0.02 : 0.15
%
% Outputs:
%   beta_map
%   Pf_FORM_map
%   Pf_MC_map
%   alphaE_map
%   alphaF_map
%   error_map
%   dominance_map
%
% ============================================================

clear; clc;
addpath(genpath(pwd));

%% ------------------------------------------------------------
% 1. Deterministic Beam Parameters
%% ------------------------------------------------------------

L = 2.0;       % [m]
b = 0.05;      % [m]
h = 0.10;      % [m]

delta_allow = 0.0032;

I = section_inertia_rect(b,h);

%% ------------------------------------------------------------
% 2. Mean Stochastic Parameters
%% ------------------------------------------------------------

mu_E = 210e9;
mu_F = 1000;

%% ------------------------------------------------------------
% 3. CoV Grid Definition
%% ------------------------------------------------------------

CoV_E_vals = 0.02 : 0.02 : 0.15;
CoV_F_vals = 0.02 : 0.02 : 0.15;

nE = length(CoV_E_vals);
nF = length(CoV_F_vals);

%% ------------------------------------------------------------
% 4. Result Matrices
%% ------------------------------------------------------------

beta_map     = zeros(nE,nF);
Pf_FORM_map  = zeros(nE,nF);
Pf_MC_map    = zeros(nE,nF);

alphaE_map   = zeros(nE,nF);
alphaF_map   = zeros(nE,nF);

error_map    = zeros(nE,nF);
dominance_map = zeros(nE,nF);

%% ------------------------------------------------------------
% 5. Monte Carlo Setup
%% ------------------------------------------------------------

N = 1e6;

%% ------------------------------------------------------------
% 6. Reliability Grid Loop
%% ------------------------------------------------------------

for i = 1:nE
for j = 1:nF

    CoV_E = CoV_E_vals(i);
    CoV_F = CoV_F_vals(j);

    sigma_E = CoV_E * mu_E;
    sigma_F = CoV_F * mu_F;

    % Parameter struct
    params.L = L;
    params.I = I;
    params.delta_allow = delta_allow;

    params.mu_E = mu_E;
    params.sigma_E = sigma_E;

    params.mu_F = mu_F;
    params.sigma_F = sigma_F;

    params.F = mu_F;

    %% ------------------------------
    % FORM (HL-RF)
    %% ------------------------------

    [beta_FORM, Z_star, Pf_FORM] = form_2d_hlrf(params);

    alpha = importance_factors(Z_star,beta_FORM);

    %% ------------------------------
    % Monte Carlo
    %% ------------------------------

    [E_samples,F_samples] = generate_EF_samples( ...
        mu_E, sigma_E, mu_F, sigma_F, N);

    [Pf_MC,~] = failure_probability_mc_2d( ...
        E_samples, F_samples, params);

    %% ------------------------------
    % Store Results
    %% ------------------------------

    beta_map(i,j) = beta_FORM;

    Pf_FORM_map(i,j) = Pf_FORM;
    Pf_MC_map(i,j)   = Pf_MC;

    alphaE_map(i,j) = alpha(1);
    alphaF_map(i,j) = alpha(2);

    error_map(i,j) = abs(Pf_MC - Pf_FORM);

    %% Dominance classification

    if abs(alpha(1)) > abs(alpha(2))
        dominance_map(i,j) = -1;   % material dominated
    else
        dominance_map(i,j) = +1;   % load dominated
    end

end
end

%% ------------------------------------------------------------
% 7. Save Results
%% ------------------------------------------------------------

results.beta      = beta_map;
results.Pf_FORM   = Pf_FORM_map;
results.Pf_MC     = Pf_MC_map;

results.alphaE    = alphaE_map;
results.alphaF    = alphaF_map;

results.error     = error_map;

results.CoV_E     = CoV_E_vals;
results.CoV_F     = CoV_F_vals;

save results_phase2C.mat results

%% ------------------------------------------------------------
% 8. Plotting
%% ------------------------------------------------------------

[CoV_F_grid,CoV_E_grid] = meshgrid(CoV_F_vals,CoV_E_vals);

%% 1 — Reliability Surface

figure
surf(CoV_F_grid,CoV_E_grid,beta_map)

xlabel('CoV_F')
ylabel('CoV_E')
zlabel('\beta')

title('Reliability Index Surface')

colorbar
shading interp
view(45,30)
grid on

%% 2 — Importance Factor Surface

figure
surf(CoV_F_grid,CoV_E_grid,alphaE_map)

xlabel('CoV_F')
ylabel('CoV_E')
zlabel('\alpha_E')

title('Importance Factor Surface (\alpha_E)')

colorbar
shading interp
view(45,30)
grid on

%% 3 — Dominance Map

figure
imagesc(CoV_F_vals,CoV_E_vals,dominance_map)

set(gca,'YDir','normal')

xlabel('CoV_F')
ylabel('CoV_E')

title('Failure Dominance Map')

colorbar
colormap(jet)

%% 4 — FORM Error Map

figure
surf(CoV_F_grid,CoV_E_grid,error_map)

xlabel('CoV_F')
ylabel('CoV_E')
zlabel('|Pf_{MC} - Pf_{FORM}|')

title('FORM Approximation Error')

colorbar
shading interp
view(45,30)
grid on