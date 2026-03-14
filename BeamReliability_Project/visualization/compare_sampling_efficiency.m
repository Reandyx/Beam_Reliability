function compare_sampling_efficiency(params)

rng(1)

% sample sizes
sample_sizes = round(logspace(3,5,6));

Pf_MC  = zeros(size(sample_sizes));
Pf_IS  = zeros(size(sample_sizes));
Pf_AIS = zeros(size(sample_sizes));

% FORM reference solution
[~,Z_star,Pf_FORM] = form_2d_hlrf(params);

for i = 1:length(sample_sizes)

    N = sample_sizes(i);

    % ---------------------------
    % Monte Carlo
    % ---------------------------

    [E,F] = generate_EF_samples( ...
        params.mu_E,params.sigma_E, ...
        params.mu_F,params.sigma_F,N);

    [Pf_MC(i),~] = failure_probability_mc_2d(E,F,params);

    % ---------------------------
    % Importance Sampling
    % ---------------------------

    Sigma = eye(2);

    [Pf_IS(i),~,~,~] = importance_sampling_pf(params,Z_star,Sigma,N);

    % ---------------------------
    % Adaptive Importance Sampling
    % ---------------------------

    max_iter = 15;

    [Pf_AIS(i),~,~] = adaptive_importance_sampling(params,Z_star,Sigma,N,max_iter);

end

% ---------------------------
% Variance estimate
% ---------------------------

var_MC = Pf_MC .* (1 - Pf_MC) ./ sample_sizes;

var_IS  = abs(Pf_IS - Pf_FORM).^2;
var_AIS = abs(Pf_AIS - Pf_FORM).^2;

% Variance reduction factor
VRF_IS  = var_MC ./ var_IS;
VRF_AIS = var_MC ./ var_AIS;

figure
hold on

loglog(sample_sizes,VRF_IS,'s-','LineWidth',2,'MarkerSize',7)
loglog(sample_sizes,VRF_AIS,'d-','LineWidth',2,'MarkerSize',7)

ylim([1 100])
xlim([min(sample_sizes) max(sample_sizes)])

xlabel('Number of samples N')
ylabel('Variance reduction factor')
title('Efficiency of importance sampling')

legend('Importance sampling','Adaptive importance sampling','Location','northwest')

grid on
box on
set(gca,'FontSize',12,'LineWidth',1.2)

end