function method_comparison_plot(params)

N = 1e5;

%% Analytical
[beta_a,Pf_a] = analytical_reliability_estimate(params);

%% FORM
[beta_FORM,~,Pf_FORM] = form_2d_hlrf(params);

%% SORM
Pf_SORM = sorm_correction(params);

%% Monte Carlo
[E_samples,F_samples] = generate_EF_samples( ...
    params.mu_E, params.sigma_E, ...
    params.mu_F, params.sigma_F, N);

[Pf_MC,~] = failure_probability_mc_2d(E_samples, F_samples, params);

%% Data
methods = {'MC','FORM','SORM','Analytical'};
values  = [Pf_MC Pf_FORM Pf_SORM Pf_a];

figure
b = bar(values,'FaceColor',[0.2 0.5 0.8]);

set(gca,'XTickLabel',methods)

ylabel('Failure Probability')
title('Reliability Method Comparison')

grid on
box on
set(gca,'FontSize',12)

%% Zoom Y-axis so differences are visible
ymin = min(values)*0.995;
ymax = max(values)*1.005;
ylim([ymin ymax])

Pf_MC_ref = values(1);

for i = 1:length(values)

    val = values(i);

    if i == 1
        label = sprintf('%.4f\n(reference)',val);
    else
        err = 100*(val - Pf_MC_ref)/Pf_MC_ref;
        label = sprintf('%.4f\n(%+.1f%% vs MC)',val,err);
    end

    text(i,val+0.0003,label,...
        'HorizontalAlignment','center',...
        'FontSize',10)

end

end