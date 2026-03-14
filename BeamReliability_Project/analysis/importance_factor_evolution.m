function importance_factor_evolution(results)

alphaE = results.alphaE;
alphaF = results.alphaF;

CoV_E = results.CoV_E;
CoV_F = results.CoV_F;

[CoV_F_grid,CoV_E_grid] = meshgrid(CoV_F,CoV_E);

%% Surface α_E

figure
surf(CoV_F_grid,CoV_E_grid,alphaE)

xlabel('CoV_F')
ylabel('CoV_E')
zlabel('\alpha_E')

title('Importance Factor Surface (\alpha_E)')
colorbar
shading interp
view(45,30)
grid on


%% Surface α_F

figure
surf(CoV_F_grid,CoV_E_grid,alphaF)

xlabel('CoV_F')
ylabel('CoV_E')
zlabel('\alpha_F')

title('Importance Factor Surface (\alpha_F)')
colorbar
shading interp
view(45,30)
grid on


%% Transition curve (CoV_E = CoV_F)

n = length(CoV_E);

alphaE_diag = zeros(n,1);
alphaF_diag = zeros(n,1);

for k = 1:n
    alphaE_diag(k) = abs(alphaE(k,k));
    alphaF_diag(k) = abs(alphaF(k,k));
end

figure

plot(CoV_E,alphaE_diag,'o-','LineWidth',2)
hold on
plot(CoV_E,alphaF_diag,'s-','LineWidth',2)

xlabel('CoV_E = CoV_F')
ylabel('|Importance Factor|')

legend('|α_E|','|α_F|')

title('Transition Between Material and Load Dominance')

grid on

%% Transition boundary |alphaE| = |alphaF|

diff_map = abs(alphaE) - abs(alphaF);

figure

[CoV_F_grid,CoV_E_grid] = meshgrid(CoV_F,CoV_E);

contour(CoV_F_grid,CoV_E_grid,diff_map,[0 0],'LineWidth',3)

xlabel('CoV_F')
ylabel('CoV_E')

title('Transition Boundary: |α_E| = |α_F|')

grid on