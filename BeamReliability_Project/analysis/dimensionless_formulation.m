function dimensionless_formulation(params)

mu_E = params.mu_E;
mu_F = params.mu_F;

L = params.L;
I = params.I;
delta_allow = params.delta_allow;

%% Dimensionless constant

C = (mu_F * L^3) / (3 * mu_E * I * delta_allow);

fprintf('\n===== Dimensionless Limit State =====\n')

fprintf('g* = 1 - C * (F_tilde / E_tilde)\n')

fprintf('\nDimensionless constant C = %.6f\n',C)

fprintf('\nNormalized variables:\n')
fprintf('E_tilde = E / mu_E\n')
fprintf('F_tilde = F / mu_F\n')

fprintf('\nGoverning parameters:\n')
fprintf('C (structural scaling)\n')
fprintf('CoV_E (material uncertainty)\n')
fprintf('CoV_F (load uncertainty)\n')

fprintf('\nInterpretation:\n')

fprintf(['Reliability depends on the ratio of normalized load ' ...
         'to normalized stiffness.\n'])

fprintf(['Increasing load uncertainty drives failure when ' ...
         '|α_F| dominates.\n'])

fprintf(['Increasing material uncertainty drives failure ' ...
         'when |α_E| dominates.\n\n'])

E_tilde = linspace(0.7,1.3,50);
F_tilde = linspace(0.7,1.3,50);

[F_grid,E_grid] = meshgrid(F_tilde,E_tilde);

g_star = 1 - C .* (F_grid ./ E_grid);

figure
surf(F_grid,E_grid,g_star)

xlabel('F̃')
ylabel('Ẽ')
zlabel('g*')

title('Dimensionless Limit State Surface')

colorbar
shading interp