function plot_mc_form_comparison(params)

N_values = logspace(2,6,15);

Pf_mc = zeros(size(N_values));

for i = 1:length(N_values)

    N = round(N_values(i));

    [Pf_mc(i),~] = failure_probability_mc_2d(params,N);

end

[beta_form,~,Pf_form] = form_2d_hlrf(params);

figure
semilogx(N_values,Pf_mc,'o-','LineWidth',2)
hold on

yline(Pf_form,'r--','LineWidth',2)

xlabel('Number of Samples')
ylabel('Failure Probability')

legend('Monte Carlo','FORM')

title('Monte Carlo Convergence vs FORM Estimate')

grid on

end