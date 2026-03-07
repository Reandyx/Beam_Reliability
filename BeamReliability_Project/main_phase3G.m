clc
clear
close all

disp('==============================')
disp('PHASE 3G — RELIABILITY DESIGN')
disp('==============================')

%% PARAMETERS

params.mu_E = 210e9;
params.sigma_E = 21e9;

params.mu_F = 1000;
params.sigma_F = 200;

params.L = 2;
params.I = 8e-6;

params.delta_allow = 0.005;

params.verbose = false;

beta_target = 3;

% Beam width
b = 0.05;

%% DESIGN CURVE

h_vals = linspace(0.02,0.15,40);

beta_vals = zeros(size(h_vals));

disp('Computing design curve...')

for i=1:length(h_vals)

    h = h_vals(i);

    I = b*h^3/12;

    params.I = I;

    [beta,~,~] = form_2d_hlrf(params);

    beta_vals(i) = beta;

end

%% SOLVE DESIGN PROBLEM

[h_opt, beta_opt, I_opt] = reliability_based_design(params,beta_target,b);

%% PRINT RESULTS

disp(' ')
disp(['Target reliability beta = ',num2str(beta_target)])

disp(' ')
disp('Optimal beam design:')

fprintf('h_opt  = %.4f m\n',h_opt)
fprintf('I_opt  = %.4e m^4\n',I_opt)
fprintf('beta   = %.4f\n',beta_opt)

%% VISUALIZATION

plot_design_curve(h_vals,beta_vals,h_opt,beta_opt,beta_target)