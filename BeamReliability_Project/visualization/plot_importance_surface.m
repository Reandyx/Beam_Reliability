function plot_importance_surface(params)

cov_E = linspace(0.01,0.3,25);
cov_F = linspace(0.01,0.3,25);

alpha_E = zeros(length(cov_E),length(cov_F));
alpha_F = zeros(length(cov_E),length(cov_F));

for i = 1:length(cov_E)

for j = 1:length(cov_F)

params.E_cov = cov_E(i);
params.F_cov = cov_F(j);

[~,U_star,~] = form_2d_hlrf(params);

alpha = importance_factors(U_star);

alpha_E(i,j) = alpha(1);
alpha_F(i,j) = alpha(2);

end
end

figure

surf(cov_E,cov_F,alpha_F')

xlabel('CoV_E')
ylabel('CoV_F')
zlabel('\alpha_F')

title('Importance Factor Surface (Load)')

shading interp
colorbar

end