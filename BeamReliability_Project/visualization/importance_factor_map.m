function importance_factor_map(params)

CoV_E = linspace(0.01,0.4,40);
CoV_F = linspace(0.01,0.4,40);

alphaE = zeros(length(CoV_E),length(CoV_F));

for i = 1:length(CoV_E)
for j = 1:length(CoV_F)

params.sigma_E = params.mu_E * CoV_E(i);
params.sigma_F = params.mu_F * CoV_F(j);

[beta,Z_star,~] = form_2d_hlrf(params);

alpha = Z_star / beta;

alphaE(i,j) = abs(alpha(1));

end
end

[CoV_F_grid,CoV_E_grid] = meshgrid(CoV_F,CoV_E);

figure
contourf(CoV_F_grid,CoV_E_grid,alphaE,20)

colorbar

xlabel('CoV_F')
ylabel('CoV_E')

title('|α_E| Importance Factor Map')

end