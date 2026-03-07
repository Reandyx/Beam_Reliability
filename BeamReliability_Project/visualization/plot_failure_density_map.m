function plot_failure_density_map(params)

% PLOT_FAILURE_DENSITY_MAP
% Visualizes probability density weighted by failure region

u_range = linspace(-6,6,200);
[U1,U2] = meshgrid(u_range,u_range);

density = zeros(size(U1));

for i = 1:numel(U1)

    Z = [U1(i); U2(i)];

    g = limit_state_EF(Z,params);

    % Standard normal density
    p = mvnpdf([U1(i) U2(i)],[0 0],eye(2));

    if g < 0
        density(i) = p;
    else
        density(i) = 0;
    end

end

figure

contourf(U1,U2,density,30,'LineColor','none')
colorbar

xlabel('U_E')
ylabel('U_F')

title('Failure Probability Density Map')

grid on

end