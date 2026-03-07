function plot_failure_region(params, Z_star, U_samples, g_vals)

% PLOT_FAILURE_REGION
% Visualizes failure region and importance sampling cloud in U-space

% Create grid in U-space
u_range = linspace(-6,6,200);
[U1,U2] = meshgrid(u_range,u_range);

g_grid = zeros(size(U1));

for i = 1:numel(U1)

    Z = [U1(i); U2(i)];

    g_grid(i) = limit_state_EF(Z,params);

end

figure

% -------------------------------------------------------
% Failure boundary (g(U)=0)
% -------------------------------------------------------
contour(U1,U2,g_grid,[0 0],'k','LineWidth',2)
hold on

% -------------------------------------------------------
% Reliability index circle (||U|| = beta)
% -------------------------------------------------------
beta = norm(Z_star);

theta = linspace(0,2*pi,300);

U_circle_x = beta * cos(theta);
U_circle_y = beta * sin(theta);

plot(U_circle_x, U_circle_y,'c--','LineWidth',1.8)

% -------------------------------------------------------
% Safe / failure samples
% -------------------------------------------------------
fail_idx = g_vals < 0;

scatter(U_samples(~fail_idx,1),U_samples(~fail_idx,2),10,'b','filled')
scatter(U_samples(fail_idx,1),U_samples(fail_idx,2),10,'r','filled')

% -------------------------------------------------------
% FORM design point
% -------------------------------------------------------
plot(Z_star(1),Z_star(2),'kp','MarkerSize',14,'LineWidth',2)

% -------------------------------------------------------
% Sampling distribution ellipse
% -------------------------------------------------------
theta = linspace(0,2*pi,200);
circle = [cos(theta); sin(theta)];

Sigma = cov(U_samples);
ellipse = chol(Sigma) * circle;

plot(Z_star(1) + ellipse(1,:), ...
     Z_star(2) + ellipse(2,:), ...
     'm','LineWidth',2)

% -------------------------------------------------------
% Labels and legend
% -------------------------------------------------------
xlabel('U_E')
ylabel('U_F')

title('Failure Region and Importance Sampling Cloud')

legend('Failure boundary','Reliability index \beta','Safe samples',...
       'Failure samples','FORM design point','Sampling distribution')

xlim([-7 -3])
ylim([1 6])
axis equal
grid on

end