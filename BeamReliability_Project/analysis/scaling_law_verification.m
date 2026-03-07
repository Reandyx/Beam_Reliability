function scaling_law_verification(params)

CoV_E = linspace(0.01,0.4,30);

beta_FORM = zeros(size(CoV_E));
beta_scaling = zeros(size(CoV_E));

CoV_F = params.sigma_F / params.mu_F;

for i = 1:length(CoV_E)

    params.sigma_E = params.mu_E * CoV_E(i);

    [beta,~,~] = form_2d_hlrf(params);
    beta_FORM(i) = beta;

    beta_raw = 1 ./ sqrt(CoV_E.^2 + CoV_F^2);

    % normalize scaling law to match first FORM value
    scale_factor = beta_FORM(1) / beta_raw(1);
    
    beta_scaling = scale_factor * beta_raw;

end

figure
plot(CoV_E,beta_FORM,'LineWidth',2)
hold on
plot(CoV_E,beta_scaling,'--','LineWidth',2)

xlabel('CoV_E')
ylabel('\beta')

legend('FORM','Scaling Law (normalized)')
title('Reliability Scaling Law Verification')

grid on

end