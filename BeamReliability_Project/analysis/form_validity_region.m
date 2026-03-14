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

    error_map(i,j) = abs(Pf_MC - Pf_FORM) / max(Pf_MC,1e-12);

end
end

%% ============================================================
% Plot FORM validity map (journal style)
%% ============================================================

% Smooth Monte Carlo noise
error_map = imgaussfilt(error_map,1.5);

[CoV_F_grid,CoV_E_grid] = meshgrid(CoV_F_vals,CoV_E_vals);

figure
hold on

%% Heatmap (fewer levels -> cleaner)
contourf(CoV_F_grid,CoV_E_grid,error_map,12,'LineColor','none')

colormap(parula)

%% Contours at meaningful error levels
levels = [0.0005 0.001 0.0015 0.002 0.0025];

[C,hc] = contour(CoV_F_grid,CoV_E_grid,error_map,levels,...
                 'k','LineWidth',1.2);

clabel(C,hc,'FontSize',10,'Color','white')

%% Colorbar
cb = colorbar;
cb.Label.String = 'Relative error (×10^{-3})';

%% Axes
xlabel('Load uncertainty CoV_F')
ylabel('Material uncertainty CoV_E')

title('FORM validity map (relative error vs Monte Carlo)')

grid on
box on
set(gca,'FontSize',12)

%% Better color scaling
caxis([0 max(error_map(:))])

%% Save results

results_dir = fullfile(pwd,'results');

if ~exist(results_dir,'dir')
    mkdir(results_dir)
end

save(fullfile(results_dir,'form_validity_error_map.mat'),'error_map')

exportgraphics(gcf,fullfile(results_dir,'phase3B_FORM_validity_map.png'),'Resolution',300)