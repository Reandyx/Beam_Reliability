function method_comparison_plot(params)

N = 1e5;

% Analytical
[beta_a,Pf_a] = analytical_reliability_estimate(params);

% FORM
[beta_FORM,~,Pf_FORM] = form_2d_hlrf(params);

% SORM
Pf_SORM = sorm_correction(params);

% Monte Carlo
[E_samples,F_samples] = generate_EF_samples( ...
    params.mu_E, params.sigma_E, ...
    params.mu_F, params.sigma_F, N);

[Pf_MC,~] = failure_probability_mc_2d(E_samples, F_samples, params);

methods = {'MC','FORM','SORM','Analytical'};
values = [Pf_MC Pf_FORM Pf_SORM Pf_a];

figure
bar(values)

set(gca,'XTickLabel',methods)

ylabel('Failure Probability')
title('Reliability Method Comparison')

grid on

end