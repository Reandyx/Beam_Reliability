function reliability_regime_diagram(params)

CoV_E = linspace(0.01,0.4,40);
CoV_F = linspace(0.01,0.4,40);

nE = length(CoV_E);
nF = length(CoV_F);

beta_map = zeros(nE,nF);
alpha_diff = zeros(nE,nF);

for i = 1:nE
for j = 1:nF

params.sigma_E = params.mu_E * CoV_E(i);
params.sigma_F = params.mu_F * CoV_F(j);

[beta,Z_star,~] = form_2d_hlrf(params);

beta_map(i,j) = beta;
Pf_map = normcdf(-beta_map);

alpha = Z_star / beta;

alpha_diff(i,j) = abs(alpha(1)) - abs(alpha(2));

end
end

[CoV_F_grid,CoV_E_grid] = meshgrid(CoV_F,CoV_E);

figure

contourf(CoV_F_grid,CoV_E_grid,beta_map,20)
colorbar
hold on

% reliability contours
contour(CoV_F_grid,CoV_E_grid,Pf_map,[0.1 0.05 0.01],'k','LineWidth',2)

% dominance boundary
contour(CoV_F_grid,CoV_E_grid,alpha_diff,[0 0],'w','LineWidth',2)

% load FORM validity map
load results/form_validity_error_map.mat

% recreate Phase 3B grid
CoV_E_vals = 0.02:0.02:0.30;
CoV_F_vals = 0.02:0.02:0.30;

[CoV_F_valid,CoV_E_valid] = meshgrid(CoV_F_vals,CoV_E_vals);

contour(CoV_F_valid,CoV_E_valid,error_map,[0.05 0.05],'r','LineWidth',3)
xlabel('CoV_F')
ylabel('CoV_E')

title('Reliability Regime Diagram')

grid on

end