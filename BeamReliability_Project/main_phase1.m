%% ============================================================
%  PHASE 1 — STOCHASTIC MODELING OF YOUNG'S MODULUS
%  Project: Probabilistic Reliability of Beam
%  Author: SA-01 (Stochastic Architect)
%  ============================================================

clear; clc;

%% ------------------ PATH SETUP ------------------
addpath(genpath(pwd));   % ensures mechanics + stochastic are visible

%% ------------------ PARAMETERS (SI UNITS ONLY) ------------------

% Young's modulus statistics
mu    = 210e9;      % Pa
sigma = 15e9;       % Pa

% Monte Carlo size
N = 1e5;

% Beam parameters
F = 5000;           % N  (adjusted to reach meaningful failure regime)
L = 2;              % m
I = 8e-6;           % m^4

% Serviceability limit
delta_allow = L / 250;   % m

%% ------------------ MONTE CARLO SIMULATION ------------------

E_samples     = generate_E_samples(mu,sigma,N);
delta_samples = monte_carlo_deflection(F,L,E_samples,I);
[Pf,CI]       = failure_probability_mc(delta_samples,delta_allow);

%% ------------------ STATISTICAL VALIDATION ------------------

mean_E = mean(E_samples);
std_E  = std(E_samples);

rel_mean_error = abs(mean_E - mu)/mu;
rel_std_error  = abs(std_E  - sigma)/sigma;

neg_rate = sum(E_samples < 0)/N;

mean_delta = mean(delta_samples);
std_delta  = std(delta_samples);

%% ------------------ PRINT PHASE 1 REPORT ------------------

fprintf('\n====================================================\n');
fprintf('PHASE 1 COMPLETION REPORT\n');
fprintf('====================================================\n');

fprintf('\n--- Young''s Modulus Statistics ---\n');
fprintf('Target mu     = %.3e Pa\n', mu);
fprintf('Sample mean   = %.3e Pa\n', mean_E);
fprintf('Target sigma  = %.3e Pa\n', sigma);
fprintf('Sample std    = %.3e Pa\n', std_E);
fprintf('Relative mean error = %.4f %%\n', rel_mean_error*100);
fprintf('Relative std error  = %.4f %%\n', rel_std_error*100);
fprintf('Negative E rate     = %.6f %%\n', neg_rate*100);

fprintf('\n--- Deflection Statistics ---\n');
fprintf('Mean(delta) = %.6e m\n', mean_delta);
fprintf('Std(delta)  = %.6e m\n', std_delta);
fprintf('Allowable   = %.6e m\n', delta_allow);

fprintf('\n--- Reliability Results ---\n');
fprintf('Failure Probability Pf = %.6f\n', Pf);
fprintf('95%% Confidence Interval = ± %.6f\n', CI);

fprintf('\n====================================================\n');

%% ------------------ REQUIRED PLOTS ------------------

% 1️ E distribution
figure
histogram(E_samples,100,'Normalization','pdf')
hold on
x_vals = linspace(min(E_samples),max(E_samples),1000);
plot(x_vals, normpdf(x_vals,mu,sigma),'LineWidth',2)
xlabel('Young''s Modulus E [Pa]')
ylabel('PDF')
title('Distribution of Young''s Modulus')
grid on
saveas(gcf,'plots/E_distribution.png')

% 2️ Deflection distribution
figure
histogram(delta_samples,100,'Normalization','pdf')
xlabel('Deflection \delta [m]')
ylabel('PDF')
title('Distribution of Beam Deflection')
grid on
saveas(gcf,'plots/delta_distribution.png')

% 3️ Failure visualization
figure
histogram(delta_samples,100)
hold on
xline(delta_allow,'r','LineWidth',2)
xlabel('Deflection \delta [m]')
ylabel('Frequency')
title('Failure Threshold Visualization')
grid on
saveas(gcf,'plots/failure_visualization.png')

%% ------------------ CONVERGENCE STUDY ------------------

convergence_study(mu,sigma,F,L,I,delta_allow);