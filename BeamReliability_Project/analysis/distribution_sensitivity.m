%% ============================================================
% PHASE 3B — DISTRIBUTION SENSITIVITY
% Compare Normal / Lognormal / Weibull material models
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

%% Load uncertainty (fixed)

CoV_F = 0.05;
sigma_F = CoV_F * mu_F;

%% Material uncertainty range

CoV_E_vals = 0.02 : 0.02 : 0.30;
nE = length(CoV_E_vals);

%% Monte Carlo

N = 1e6;

Pf_normal    = zeros(nE,1);
Pf_lognormal = zeros(nE,1);
Pf_weibull   = zeros(nE,1);

%% Loop over material uncertainty

for i = 1:nE

    CoV_E = CoV_E_vals(i);
    sigma_E = CoV_E * mu_E;

    %% Load samples (same for all cases)

    F_samples = mu_F + sigma_F .* randn(N,1);
    F_samples(F_samples <= 0) = mu_F;

    %% ---------- Normal distribution ----------

    E_samples = mu_E + sigma_E .* randn(N,1);
    E_samples(E_samples <= 0) = mu_E;

    delta = beam_deflection(F_samples,L,E_samples,I);

    Pf_normal(i) = sum(delta > delta_allow) / N;

    %% ---------- Lognormal distribution ----------

    E_samples = generate_E_samples_lognormal(mu_E,sigma_E,N);

    delta = beam_deflection(F_samples,L,E_samples,I);

    Pf_lognormal(i) = sum(delta > delta_allow) / N;

    %% ---------- Weibull distribution ----------

    E_samples = generate_E_samples_weibull(mu_E,sigma_E,N);

    delta = beam_deflection(F_samples,L,E_samples,I);

    Pf_weibull(i) = sum(delta > delta_allow) / N;

end

%% Plot comparison

figure

plot(CoV_E_vals,Pf_normal,'-o','LineWidth',2)
hold on

plot(CoV_E_vals,Pf_lognormal,'-s','LineWidth',2)
plot(CoV_E_vals,Pf_weibull,'-d','LineWidth',2)

xlabel('CoV_E')
ylabel('Failure Probability')

title('Distribution Sensitivity')

legend('Normal','Lognormal','Weibull','Location','northwest')

grid on

%% Save figure

if ~exist('results','dir')
    mkdir results
end

saveas(gcf,'results/phase3B_distribution_sensitivity.png')