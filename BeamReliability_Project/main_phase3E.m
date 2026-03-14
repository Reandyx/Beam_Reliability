clear
clc

fprintf('==============================\n')
fprintf('PHASE 3E — IMPORTANCE SAMPLING\n')
fprintf('==============================\n')

%% PARAMETERS

params.mu_E = 210e9;
params.sigma_E = 21e9;

params.mu_F = 1000;
params.sigma_F = 200;

params.L = 2;
params.I = 8e-6;

params.delta_allow = 0.005;

params.verbose = false;

%% 1️⃣ FORM ANALYSIS

fprintf('\nRunning FORM...\n')

[beta, Z_star, Pf_FORM] = form_2d_hlrf(params);

fprintf('FORM beta = %.4f\n',beta)
fprintf('FORM Pf   = %.3e\n',Pf_FORM)

%% 2️⃣ SORM CURVATURE

fprintf('\nComputing SORM curvature...\n')

kappa = compute_sorm_curvature(Z_star,params);

fprintf('Curvature kappa = %.4e\n',kappa)

Sigma = sorm_guided_covariance(kappa,2);

%% 3️⃣ MONTE CARLO BASELINE

fprintf('\nRunning Monte Carlo baseline...\n')

N_MC = 1e6;

[E_samples,F_samples] = generate_EF_samples( ...
    params.mu_E,params.sigma_E,...
    params.mu_F,params.sigma_F,N_MC);

[Pf_MC,CI_MC] = failure_probability_mc_2d(E_samples,F_samples,params);

var_MC = Pf_MC*(1-Pf_MC)/N_MC;

fprintf('Monte Carlo Pf = %.3e\n',Pf_MC)

%% 4️⃣ BASELINE IMPORTANCE SAMPLING

fprintf('\nRunning Importance Sampling...\n')

N_IS = 1e4;

[Pf_IS,var_IS,ESS,U_samples,g_vals,w] = ...
    importance_sampling_pf(params,Z_star,Sigma,N_IS);

fprintf('IS Pf        = %.3e\n',Pf_IS)
fprintf('IS variance  = %.3e\n',var_IS)
fprintf('IS ESS       = %.1f\n',ESS)

%% 5️⃣ ADAPTIVE IMPORTANCE SAMPLING

fprintf('\nRunning Adaptive Importance Sampling...\n')

AIS_iter = 5;

[Pf_AIS,var_AIS,mu_history,ESS_AIS] = ...
    adaptive_importance_sampling(params,Z_star,Sigma,N_IS,AIS_iter);

fprintf('AIS Pf       = %.3e\n',Pf_AIS)
fprintf('AIS variance = %.3e\n',var_AIS)
fprintf('AIS ESS      = %.1f\n',ESS_AIS)

%% 6️⃣ PERFORMANCE METRICS

VRF_AIS = variance_reduction_factor(var_IS,var_AIS);

fprintf('\nVariance Reduction (AIS vs IS) = %.2fx\n',VRF_AIS)
%% 7️⃣ VISUALIZATION

fprintf('\nGenerating visualizations...\n')

% Failure region geometry
plot_failure_region(params,Z_star,U_samples,g_vals)

% Failure probability density
plot_failure_density_map(params)

% AIS trajectory
figure
plot(mu_history(:,1),mu_history(:,2),'o-','LineWidth',2)

xlabel('U_E')
ylabel('U_F')

title('Adaptive Importance Sampling Trajectory')

grid on

%% DONE

fprintf('\nPhase 3E completed successfully.\n')