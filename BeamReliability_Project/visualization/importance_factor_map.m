function importance_factor_map(params)

%% Uncertainty grid
CoV_E_vals = 0.02:0.02:0.40;
CoV_F_vals = 0.02:0.02:0.40;

nE = length(CoV_E_vals);
nF = length(CoV_F_vals);

alpha_E_map = zeros(nE,nF);

%% Compute importance factors

for i = 1:nE
for j = 1:nF

    CoV_E = CoV_E_vals(i);
    CoV_F = CoV_F_vals(j);

    sigma_E = CoV_E * params.mu_E;
    sigma_F = CoV_F * params.mu_F;

    params.sigma_E = sigma_E;
    params.sigma_F = sigma_F;

    [beta,Z_star,~] = form_2d_hlrf(params);

    alpha_E_map(i,j) = abs(Z_star(1))/norm(Z_star);

end
end

%% Grid
[CoV_F_grid,CoV_E_grid] = meshgrid(CoV_F_vals,CoV_E_vals);

figure
hold on

%% Heatmap
contourf(CoV_F_grid,CoV_E_grid,alpha_E_map,80,'LineColor','none')
colormap(parula)

%% Contours
levels = 0.1:0.1:0.9;

[C,hc] = contour(CoV_F_grid,CoV_E_grid,alpha_E_map,levels,'k','LineWidth',1);
clabel(C,hc,'FontSize',10)

%% Equal influence line
contour(CoV_F_grid,CoV_E_grid,alpha_E_map,[0.707 0.707],...
        'w--','LineWidth',2)

text(0.28,0.33,'|\alpha_E| = |\alpha_F|','Color','white','FontSize',11)

%% Reference point
CoV_F_ref = 0.20;
CoV_E_ref = 0.10;

plot(CoV_F_ref,CoV_E_ref,'ko','MarkerSize',8,'MarkerFaceColor','white')
text(CoV_F_ref+0.01,CoV_E_ref,'Reference case','FontSize',10)

%% Colorbar
cb = colorbar;
cb.Label.String = '|\alpha_E| importance factor';

caxis([0 1])

%% Labels
xlabel('CoV_F')
ylabel('CoV_E')

title('|\alpha_E| Importance Factor Map')

grid on
box on
set(gca,'FontSize',12)

end