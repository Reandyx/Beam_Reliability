%% ============================================================
% GENERATE ALL PAPER FIGURES
% Probabilistic Reliability Analysis of a Cantilever Beam
% ============================================================

clear; clc; close all;

addpath(genpath(pwd));

rng(1)   % reproducibility

%% ------------------------------------------------------------
% Force publication (white) figure style
%% ------------------------------------------------------------
% 
% set(groot,'defaultFigureColor','w')
% set(groot,'defaultAxesColor','w')
% set(groot,'defaultAxesXColor','k')
% set(groot,'defaultAxesYColor','k')
% set(groot,'defaultAxesGridColor',[0.8 0.8 0.8])
% set(groot,'defaultAxesMinorGridColor',[0.9 0.9 0.9])
% set(groot,'defaultAxesFontSize',12)
% set(groot,'defaultLineLineWidth',1.8)

%% ------------------------------------------------------------
% Create output folder
%% ------------------------------------------------------------

if ~exist('paper_figures','dir')
    mkdir paper_figures
end

fprintf('\nGenerating all paper figures...\n\n')

%% ------------------------------------------------------------
% Beam parameters
%% ------------------------------------------------------------

params.L = 2.0;
params.b = 0.05;
params.h = 0.10;

params.delta_allow = 0.004; 

params.I = section_inertia_rect(params.b,params.h);

params.mu_E = 210e9;
params.mu_F = 1000;

params.sigma_E = 0.10 * params.mu_E;
params.sigma_F = 0.20* params.mu_F;

params.F = params.mu_F;

%% ============================================================
% FIGURE 3 — FORM Reliability Geometry
%% ============================================================

fprintf('Figure 3: FORM reliability geometry\n')

reliability_geometry_diagram(params)

exportgraphics(gcf,'paper_figures/Fig3_FORM_geometry.png','Resolution',300)

%% ============================================================
% FIGURE 4 — Sobol Sensitivity
%% ============================================================

fprintf('Figure 4: Sobol sensitivity\n')

N = 1e5;

[S_E,S_F,S_E_total,S_F_total] = sobol_indices(params,N);

plot_sobol_results(S_E,S_F,S_E_total,S_F_total)

exportgraphics(gcf,'paper_figures/Fig4_Sobol.png','Resolution',300)

%% ============================================================
% FIGURE 5 — Reliability Atlas
%% ============================================================

fprintf('Figure 5: Reliability atlas\n')

reliability_atlas(params)

exportgraphics(gcf,'paper_figures/Fig5_ReliabilityAtlas.png','Resolution',300)


%% ============================================================
% FIGURE 6 — FORM Validity Map
%% ============================================================

fprintf('Figure 6: FORM validity region\n')

run('analysis/form_validity_region.m')

exportgraphics(gcf,'paper_figures/Fig6_FORM_validity.png','Resolution',300)


%% ============================================================
% FIGURE 7 — Importance Factor Map
%% ============================================================

fprintf('Figure 7: Importance factor map\n')
importance_factor_map(params)

saveas(gcf,'paper_figures/Fig7_ImportanceMap.png')
exportgraphics(gcf,'paper_figures/Fig7_ImportanceMap.png','Resolution',300)


%% ============================================================
% FIGURE 8 — Reliability Regime Diagram
%% ============================================================

fprintf('Figure 8: Reliability regime diagram\n')

reliability_regime_diagram(params)

saveas(gcf,'paper_figures/Fig8_RegimeDiagram.png')
exportgraphics(gcf,'paper_figures/Fig8_RegimeDiagram.png','Resolution',300)


%% ============================================================
% FIGURE 9 — Sampling Efficiency
%% ============================================================

fprintf('Figure 9: Sampling efficiency\n')

compare_sampling_efficiency(params)

exportgraphics(gcf,'paper_figures/Fig9_SamplingEfficiency.png','Resolution',300)

%% ============================================================
% FIGURE 10 — Method Comparison
%% ============================================================

fprintf('Figure 10: Method comparison\n')

method_comparison_plot(params)

exportgraphics(gcf,'paper_figures/Fig10_MethodComparison.png','Resolution',300)


%% ============================================================
% FIGURE 11 — Reliability-Based Design Curve
%% ============================================================

fprintf('Figure 11: Reliability-based design\n')

beta_target = 3;
b = 0.05;

% Beam height sweep
h_vals = linspace(0.03,0.25,300);

beta_vals = NaN(size(h_vals));
g_vals    = NaN(size(h_vals));

% Ensure consistent stochastic parameters
params.sigma_F = 0.20 * params.mu_F;
params.sigma_E = 0.10 * params.mu_E;

for i = 1:length(h_vals)

    h = h_vals(i);

    params.h = h;
    params.I = section_inertia_rect(b,h);

    try
        [beta,~,~] = form_2d_hlrf(params);
        beta_vals(i) = beta;
    catch
        beta_vals(i) = NaN;
    end

    g_vals(i) = limit_state_EF([0;0],params);

end

%% Deterministic limit

idx_det = find(g_vals > 0,1);
h_det = h_vals(idx_det);

fprintf('Deterministic limit: h ≈ %.4f m\n',h_det)

%% Optimal design

[h_opt,beta_opt,I_opt] = reliability_based_design(params,beta_target,b);

fprintf('Optimal design: h = %.4f m   beta = %.3f\n',h_opt,beta_opt)

%% Filter valid region

valid = (~isnan(beta_vals)) & (g_vals > 0);

h_vals = h_vals(valid);
beta_vals = beta_vals(valid);

%% Plot

figure
hold on

ymax = 3.5;

% Unsafe region
patch([min(h_vals) max(h_vals) max(h_vals) min(h_vals)],...
      [0 0 beta_target beta_target],...
      [1 0.9 0.9],'EdgeColor','none','FaceAlpha',0.35)

% Safe region
patch([min(h_vals) max(h_vals) max(h_vals) min(h_vals)],...
      [beta_target beta_target ymax ymax],...
      [0.9 1 0.9],'EdgeColor','none','FaceAlpha',0.35)

% Deterministic failure region
patch([min(h_vals) h_det h_det min(h_vals)],...
      [0 0 ymax ymax],...
      [0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.6)

% Reliability curve
h_curve = plot(h_vals,beta_vals,...
               'Color',[0.85 0.55 0],...
               'LineWidth',2.8);

% Target reliability
h_target = yline(beta_target,'--r','\beta_{target}','LineWidth',1.5);

% Deterministic limit
h_detline = xline(h_det,'--k','Deterministic limit','LineWidth',1.5);

% Optimal design
h_optpt = plot(h_opt,beta_opt,'ko',...
               'MarkerSize',10,...
               'MarkerFaceColor','k');

%% Labels

text(h_det/2,0.35,'Deterministic failure','FontSize',11)

text(mean(h_vals),1.5,'UNSAFE (\beta < \beta_{target})',...
     'FontSize',11,'HorizontalAlignment','center')

text(mean(h_vals),3.2,'SAFE (\beta \geq \beta_{target})',...
     'FontSize',11,'HorizontalAlignment','center')

text(h_opt+0.004,beta_opt+0.05,...
    sprintf('h_{opt} = %.3f m\nI_{opt} = %.2e m^4',h_opt,I_opt),...
    'FontSize',11,'VerticalAlignment','bottom')

xlabel('Beam height h [m]')
ylabel('Reliability Index \beta')

title('Reliability-Based Design Curve')

legend([h_curve h_target h_detline h_optpt],...
{'Reliability curve \beta(h)',...
'Target reliability','Deterministic limit','Optimal design'},...
'Location','northwest')

grid on
box on
set(gca,'FontSize',12)

xlim([0.03 0.25])
ylim([0 ymax])

exportgraphics(gcf,'paper_figures/Fig11_RBD.png','Resolution',300)

%% ============================================================
% FIGURE 12 — Reliability-Based Design Map
%% ============================================================

fprintf('Figure 12: RBD design map\n')

h_vals = linspace(0.09,0.15,120);
covF_vals = linspace(0.05,0.30,120);

beta_map = NaN(length(covF_vals),length(h_vals));

params.sigma_E = 0.10 * params.mu_E;

for i = 1:length(covF_vals)

    params.sigma_F = covF_vals(i) * params.mu_F;

    for j = 1:length(h_vals)

        h = h_vals(j);

        params.h = h;
        params.I = section_inertia_rect(b,h);

        g_det = limit_state_EF([0;0],params);

        if g_det < 0
            beta_map(i,j) = -1;
            continue
        end

        try
            [beta,~,~] = form_2d_hlrf(params);
            beta = min(beta,6);
            beta_map(i,j) = beta;
        catch
            beta_map(i,j) = NaN;
        end

    end
end

beta_map = imgaussfilt(beta_map,1.0);

%% Deterministic limit

L = params.L;
F = params.mu_F;
E = params.mu_E;
b = 0.05;
delta_allow = params.delta_allow;

h_det = ((4*F*L^3)/(E*b*delta_allow))^(1/3);

%% Optimal design

params.sigma_F = 0.20 * params.mu_F;
params.sigma_E = 0.10 * params.mu_E;

beta_target = 3;

[h_opt,beta_opt,I_opt] = reliability_based_design(params,beta_target,b);

covF_opt = params.sigma_F / params.mu_F;

%% Plot

figure
hold on

contourf(h_vals,covF_vals,beta_map,80,'LineColor','none','HandleVisibility','off')

[C,hc] = contour(h_vals,covF_vals,beta_map,[1 2 3 4],...
                 'k','LineWidth',1,'HandleVisibility','off');

hLabels = clabel(C,hc,'FontSize',13,'LabelSpacing',350,'Interpreter','tex');

set(hLabels,'Color','k','FontWeight','bold')

% beta = 3 frontier
contour(h_vals,covF_vals,beta_map,[3 3],...
        'k--','LineWidth',3,'HandleVisibility','off')

% failure boundary
contour(h_vals,covF_vals,beta_map,[0 0],...
        'r','LineWidth',2,'HandleVisibility','off')
% deterministic limit
xline(h_det,'--r','LineWidth',2,'HandleVisibility','off')

yl = ylim;
y_center = mean(yl);

text(h_det+0.001,y_center,'Deterministic limit',...
     'Color','white','FontSize',11,...
     'Rotation',90,'VerticalAlignment','middle')

% optimal design
plot(h_opt,covF_opt,'ko',...
     'MarkerSize',9,'MarkerFaceColor','k', 'HandleVisibility','off')

text(h_opt+0.003,covF_opt,'Optimal design',...
'Color','white','FontSize',11,...
'FontWeight','bold',...
'BackgroundColor','black',...
'Margin',3)

%% Region labels

text(0.138,0.09,'SAFE (\beta > 3)',...
'Color','white','FontSize',12,...
'FontWeight','bold',...
'BackgroundColor','black','Margin',3)

text(0.10,0.26,'UNSAFE (\beta < 3)',...
     'Color','white','FontSize',12,'FontWeight','bold')

%% Formatting

colormap(parula)

cb = colorbar;
cb.Label.String = 'Reliability index \beta';

caxis([0 4])

xlabel('Beam height h [m]')
ylabel('Load uncertainty CoV_F')

title('Reliability-Based Design Map')

grid on
box on
set(gca,'FontSize',12)

plot(nan,nan,'k-','LineWidth',1,'DisplayName','\beta contours');
plot(nan,nan,'k--','LineWidth',3,'DisplayName','\beta = 3 frontier');
plot(nan,nan,'r-','LineWidth',2,'DisplayName','Failure boundary');
plot(nan,nan,'r--','LineWidth',2,'DisplayName','Deterministic limit');
plot(nan,nan,'ko','MarkerFaceColor','k','MarkerSize',9,...
     'DisplayName','Optimal design');


legend('Location','northwest')

exportgraphics(gcf,'paper_figures/Fig12_RBD_map.png','Resolution',300)

%% ============================================================
% FIGURE 13 — Reliability Index vs Failure Probability
%% ============================================================

fprintf('Figure: Reliability index vs failure probability\n')

beta = linspace(0,5,300);
Pf = normcdf(-beta);

figure
semilogy(beta,Pf,'LineWidth',2)
hold on

%% Target reliability
beta_target = 3;
Pf_target = normcdf(-beta_target);

plot(beta_target,Pf_target,'ro','MarkerFaceColor','r','MarkerSize',8)

xline(beta_target,'--r','Target reliability','LineWidth',1.5)
yline(Pf_target,'--r','LineWidth',1.5)

text(beta_target+0.25, Pf_target*2.2, ...
    sprintf('\\beta = %.1f   P_f = %.2e', beta_target, Pf_target), ...
    'FontSize',11, ...
    'BackgroundColor','w', ...
    'EdgeColor',[0.8 0.8 0.8], ...
    'Margin',4)

%% Standard reliability levels
beta_levels = [1 2 3 4];
Pf_levels = normcdf(-beta_levels);

plot(beta_levels,Pf_levels,'ko','MarkerFaceColor','k')

for i = 1:length(beta_levels)

    text(beta_levels(i)+0.05,Pf_levels(i)*1.2,...
        sprintf('\\beta = %d',beta_levels(i)),...
        'FontSize',10)

end

%% Axis formatting

xlabel('Reliability index \beta')
ylabel('Failure probability P_f')

title('Reliability Index vs Failure Probability','FontWeight','normal')
grid on
grid minor

set(gca,'YScale','log')

xlim([0 5])
ylim([1e-7 1])

set(gca,'FontSize',12)

box on

exportgraphics(gcf,'paper_figures/Fig13RIvsFP.png','Resolution',300)


%% ============================================================
% FIGURE 14 — FORM vs MC convergence
%% ============================================================
disp('Figure 14: FORM vs MC convergence')
disp(Pf_FORM)
disp(Pf_MC)
form_vs_mc_convergence(params)
exportgraphics(gcf,'paper_figures/Fig14FORMvsMCcov.png','Resolution',300)

%% ============================================================
% FIGURE 15 — reliability_sensitivity_covF
%% ============================================================
disp('Figure 15: reliability_sensitivity_covF')
reliability_sensitivity_covF(params)
exportgraphics(gcf,'paper_figures/Fig15ReliabilitySensitivityCovF.png','Resolution',300)
