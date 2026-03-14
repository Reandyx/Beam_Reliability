function reliability_sensitivity_covF(params)

%% Range of load uncertainty
covF = linspace(0.05,0.40,40);

beta_vals = zeros(size(covF));
Pf_vals   = zeros(size(covF));

for i = 1:length(covF)

    params.sigma_F = covF(i) * params.mu_F;

    [beta,~,Pf] = form_2d_hlrf(params);

    beta_vals(i) = beta;
    Pf_vals(i)   = Pf;

end


%% Plot
figure
hold on
set(gcf,'Color','w')

plot(covF,beta_vals,'LineWidth',2.5)

xlabel('Load uncertainty CoV_F')
ylabel('Reliability index \beta')

title('Reliability sensitivity to load uncertainty')

grid on
box on
set(gca,'FontSize',12,'LineWidth',1)

end