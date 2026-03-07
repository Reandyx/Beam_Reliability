%% ============================================================
% PHASE 3B — FORM VALIDITY REGION
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

%% Extended uncertainty grid

CoV_E_vals = 0.02 : 0.02 : 0.30;
CoV_F_vals = 0.02 : 0.02 : 0.30;

nE = length(CoV_E_vals);
nF = length(CoV_F_vals);

beta_map = zeros(nE,nF);
Pf_FORM_map = zeros(nE,nF);
Pf_MC_map = zeros(nE,nF);

error_map = zeros(nE,nF);

%% Monte Carlo

N = 1e6;

for i = 1:nE
for j = 1:nF

    CoV_E = CoV_E_vals(i);
    CoV_F = CoV_F_vals(j);

    sigma_E = CoV_E * mu_E;
    sigma_F = CoV_F * mu_F;

    params.L = L;
    params.I = I;
    params.delta_allow = delta_allow;

    params.mu_E = mu_E;
    params.sigma_E = sigma_E;

    params.mu_F = mu_F;
    params.sigma_F = sigma_F;

    params.F = mu_F;

    %% FORM

    [beta_FORM, Z_star, Pf_FORM] = form_2d_hlrf(params);

    %% Monte Carlo

    [E_samples,F_samples] = generate_EF_samples( ...
        mu_E,sigma_E,mu_F,sigma_F,N);

    [Pf_MC,~] = failure_probability_mc_2d( ...
        E_samples,F_samples,params);

    %% Store

    beta_map(i,j) = beta_FORM;

    Pf_FORM_map(i,j) = Pf_FORM;
    Pf_MC_map(i,j)   = Pf_MC;

    error_map(i,j) = abs(Pf_MC - Pf_FORM) / Pf_MC;

end
end

%% Plot FORM error surface

[CoV_F_grid,CoV_E_grid] = meshgrid(CoV_F_vals,CoV_E_vals);

figure

% 3D error surface
surf(CoV_F_grid,CoV_E_grid,error_map)

xlabel('CoV_F')
ylabel('CoV_E')
zlabel('Relative Error')

title('FORM Validity Map')

colorbar
shading interp
view(45,30)
grid on
hold on

% 5% error boundary (FORM validity limit)
contour3(CoV_F_grid,CoV_E_grid,error_map,[0.05 0.05],'r','LineWidth',2)

%% Save results

if ~exist('results','dir')
    mkdir results
end

save('results/form_validity_error_map.mat','error_map')

saveas(gcf,'results/phase3B_FORM_validity_map.png')