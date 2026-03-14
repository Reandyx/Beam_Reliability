%% ============================================================
% PHASE 3B — CORRELATED UNCERTAINTY ANALYSIS
% Study effect of correlation between E and F
% ============================================================

clear; clc;
addpath(genpath(pwd));

%% Beam parameters

L = 2.0;
b = 0.05;
h = 0.10;

delta_allow = 0.0032;

I = section_inertia_rect(b,h);

%% Means

mu_E = 210e9;
mu_F = 1000;

%% Uncertainty levels

CoV_E = 0.10;
CoV_F = 0.10;

sigma_E = CoV_E * mu_E;
sigma_F = CoV_F * mu_F;

%% Correlation values

rho_vals = [-0.5 0 0.5];

%% Monte Carlo

N = 1e6;

Pf_vals = zeros(length(rho_vals),1);

%% Loop over correlation

for k = 1:length(rho_vals)

    rho = rho_vals(k);

    %% Correlation matrix

    R = [1 rho;
         rho 1];

    Lc = chol(R,'lower');

    %% Generate independent normals

    Z = randn(N,2);

    %% Apply correlation

    Z_corr = Z * Lc';

    Z_E = Z_corr(:,1);
    Z_F = Z_corr(:,2);

    %% Transform to physical space

    E_samples = mu_E + sigma_E * Z_E;
    F_samples = mu_F + sigma_F * Z_F;

    %% Remove nonphysical values

    E_samples(E_samples <= 0) = mu_E;
    F_samples(F_samples <= 0) = mu_F;

    %% Compute deflection

    delta = beam_deflection(F_samples,L,E_samples,I);

    %% Failure probability

    Pf_vals(k) = sum(delta > delta_allow) / N;

end

%% Plot

figure

plot(rho_vals,Pf_vals,'-o','LineWidth',2)

xlabel('\rho(E,F)')
ylabel('Failure Probability')

title('Effect of Correlation Between E and F')

grid on

%% Save figure

if ~exist('results','dir')
    mkdir results
end

saveas(gcf,'results/phase3B_correlation_effect.png')