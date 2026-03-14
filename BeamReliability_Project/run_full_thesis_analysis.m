clear
clc
close all

addpath(genpath(pwd))

disp('========================================')
disp('CANTILEVER BEAM RELIABILITY STUDY')
disp('FULL THESIS PIPELINE')
disp('========================================')

figure_id = 1;

%% PARAMETERS

params.mu_E = 210e9;
params.sigma_E = 21e9;

params.mu_F = 5000;
params.sigma_F = 1000;

params.L = 2;

params.b = 0.05;
params.I = 8e-6;

params.delta_allow = 0.01;

params.verbose = false;

N = 2e5;

%% =========================================================
% FIGURE 1 — BEAM MODEL
%% =========================================================

figure(figure_id); figure_id = figure_id + 1;

L = params.L;

plot([0 L],[0 0],'k','LineWidth',4)
hold on

plot(0,0,'ks','MarkerFaceColor','k','MarkerSize',10)
plot(L,0,'ro','MarkerFaceColor','r','MarkerSize',8)

text(L/2,0.05,'Cantilever Beam','HorizontalAlignment','center')
text(L,0.08,'Load F','HorizontalAlignment','center')

axis equal
axis([-0.2 L+0.2 -0.2 0.2])

title('Cantilever Beam Model')

grid on


%% =========================================================
% MONTE CARLO BASELINE
%% =========================================================

disp('Running Monte Carlo simulation')

[E_samples,F_samples] = generate_EF_samples( ...
    params.mu_E, ...
    params.sigma_E, ...
    params.mu_F, ...
    params.sigma_F, ...
    N);

Z_E = (E_samples - params.mu_E) ./ params.sigma_E;
Z_F = (F_samples - params.mu_F) ./ params.sigma_F;

g_vals = zeros(N,1);

for i=1:N

    Z = [Z_E(i); Z_F(i)];

    g_vals(i) = limit_state_EF(Z,params);

end

Pf_MC = mean(g_vals < 0);

CI = 1.96*sqrt(Pf_MC*(1-Pf_MC)/N);

beta_MC = -norminv(Pf_MC);

fprintf('\nMonte Carlo Results\n')
fprintf('Pf = %.6e\n',Pf_MC)
fprintf('beta = %.4f\n',beta_MC)
fprintf("E range: %.3e  %.3e\n",min(E_samples),max(E_samples))
fprintf("F range: %.3e  %.3e\n",min(F_samples),max(F_samples))

%% =========================================================
% FORM
%% =========================================================

disp('Running FORM (HLRF)')

[beta_FORM,Z_star,Pf_FORM,Z_hist] = form_2d_hlrf(params);

fprintf('\nFORM Results\n')
fprintf('Pf = %.6e\n',Pf_FORM)
fprintf('beta = %.4f\n',beta_FORM)

%% =========================================================
% IMPORTANCE FACTORS
%% =========================================================

alpha = importance_factors(beta_FORM,Z_star);

fprintf('\nImportance factors\n')
disp(alpha)

%% =========================================================
% SOBOL GLOBAL SENSITIVITY
%% =========================================================

disp('Running Sobol sensitivity')

[S_E,S_F,S_E_total,S_F_total] = sobol_indices(params,50000);

figure(figure_id); figure_id = figure_id + 1;

bar([S_E S_F; S_E_total S_F_total]')

legend('First Order','Total Order')

xticklabels({'E','F'})
ylabel('Sobol Index')

title('Global Sensitivity Analysis')

grid on

%% =========================================================
% IMPORTANCE SAMPLING
%% =========================================================

disp('Running Importance Sampling')

Sigma = eye(2);

[Pf_IS,var_IS,ESS,U,g_vals_IS,w] = importance_sampling_pf(params,Z_star,Sigma,50000);

fprintf('\nImportance Sampling Pf = %.6e\n',Pf_IS)
fprintf('Variance = %.3e\n',var_IS)
fprintf('ESS = %.1f\n',ESS)

%% =========================================================
% FIGURE — METHOD COMPARISON
%% =========================================================

figure(figure_id); figure_id = figure_id + 1;

bar([Pf_MC Pf_FORM Pf_IS])

xticklabels({'Monte Carlo','FORM','Importance Sampling'})

ylabel('Failure Probability')

title('Method Comparison')

grid on

%% =========================================================
% RELIABILITY REGIME MAP
%% =========================================================

figure(figure_id); figure_id = figure_id + 1;

try

    reliability_regime_diagram(params)

    title('Reliability Regime Diagram')

catch

    warning('Reliability regime diagram failed')

end

%% =========================================================
% RELIABILITY ATLAS
%% =========================================================

figure(figure_id); figure_id = figure_id + 1;

warning('off','all')
reliability_atlas(params)
warning('on','all')

%% =========================================================
% MONTE CARLO CONVERGENCE
%% =========================================================

N_vec = logspace(3,6,15);
Pf_conv = zeros(size(N_vec));

for i=1:length(N_vec)

    Ntest = round(N_vec(i));

    [E,F] = generate_EF_samples( ...
        params.mu_E, ...
        params.sigma_E, ...
        params.mu_F, ...
        params.sigma_F, ...
        Ntest);

    Pf_conv(i) = failure_probability_mc_2d(E,F,params);

end

figure(figure_id); figure_id = figure_id + 1;

semilogx(N_vec,Pf_conv,'o-','LineWidth',2)
hold on

yline(Pf_FORM,'r--','LineWidth',2)

xlabel('Samples')
ylabel('Failure Probability')

legend('Monte Carlo','FORM')

title('Monte Carlo Convergence')

grid on

%% =========================================================
% FAILURE REGION VISUALIZATION
%% =========================================================

[E_plot,F_plot] = generate_EF_samples( ...
    params.mu_E, ...
    params.sigma_E, ...
    params.mu_F, ...
    params.sigma_F, ...
    120000);

g = zeros(length(E_plot),1);

for i=1:length(E_plot)

    Z = [
        (E_plot(i)-params.mu_E)/params.sigma_E
        (F_plot(i)-params.mu_F)/params.sigma_F
        ];

    g(i) = limit_state_EF(Z,params);

end

figure(figure_id); figure_id = figure_id + 1;

scatter(E_plot(g>0),F_plot(g>0),5,'b','filled')
hold on

scatter(E_plot(g<0),F_plot(g<0),5,'r','filled')

xlabel('Young Modulus')
ylabel('Load')

legend('Safe','Failure')

title('Failure Region in Parameter Space')

grid on

%% =========================================================
% RELIABILITY BASED DESIGN
%% =========================================================

disp('Running reliability-based design')

beta_target = 3;

[h_opt,beta_opt,I_opt] = reliability_based_design(params,beta_target,params.b);

fprintf('\nOptimal height = %.4f m\n',h_opt)
fprintf('beta achieved = %.4f\n',beta_opt)

%% DESIGN CURVE

h_vals = linspace(0.04,0.20,40);

beta_vals = zeros(size(h_vals));

for i=1:length(h_vals)

    h = h_vals(i);

    I = params.b*h^3/12;

    params.I = I;

    try

        [beta,~,~] = form_2d_hlrf(params);

        beta_vals(i) = beta;

    catch

        beta_vals(i) = NaN;

    end

end

figure(figure_id); figure_id = figure_id + 1;

plot_design_curve(h_vals,beta_vals,h_opt,beta_opt,beta_target)

%% =========================================================
% SAVE RESULTS
%% =========================================================

results.Pf_MC = Pf_MC;
results.Pf_FORM = Pf_FORM;
results.Pf_IS = Pf_IS;

results.beta_FORM = beta_FORM;
results.beta_MC = beta_MC;

results.alpha = alpha;

results.Sobol = [S_E S_F S_E_total S_F_total];

results.h_opt = h_opt;

save results_final.mat results

disp(' ')
disp('====================================')
disp('FULL STUDY COMPLETE')
disp('====================================')