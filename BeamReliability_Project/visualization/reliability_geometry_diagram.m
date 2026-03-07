function reliability_geometry_diagram(params)

fprintf('Generating reliability geometry diagram...\n')

% --- FORM result ---
[beta, Z_star, Pf] = form_2d_hlrf(params);

fprintf('beta = %.3f\n',beta)
fprintf('Pf   = %.3e\n',Pf)

% --- U-space grid ---
u = linspace(-6,6,300);
[U1,U2] = meshgrid(u,u);

G = zeros(size(U1));

for i = 1:numel(U1)

    E = params.mu_E + params.sigma_E * U1(i);
    F = params.mu_F + params.sigma_F * U2(i);

    delta = (F*params.L^3)/(3*E*params.I);

    G(i) = params.delta_allow - delta;

end

% --- Figure ---
figure
hold on

% Safe / failure shading
contourf(U1,U2,G,[min(G(:)) 0 max(G(:))],'LineStyle','none')

colormap([0.92 0.92 0.92; 0.15 0.15 0.15])

% Failure boundary
contour(U1,U2,G,[0 0],'k','LineWidth',2)

% Origin
plot(0,0,'ko','MarkerFaceColor','k','MarkerSize',7)

% Design point
plot(Z_star(1),Z_star(2),'ro','MarkerFaceColor','r','MarkerSize',8)

% Reliability vector (origin → design point)
plot([0 Z_star(1)],[0 Z_star(2)],'r','LineWidth',2)

% Importance sampling cloud
N = 700;
samples = mvnrnd(Z_star',eye(2),N);

scatter(samples(:,1),samples(:,2),10,'b','filled','MarkerFaceAlpha',0.25)

% Labels
text(Z_star(1)+0.3,Z_star(2),'Design point U*','FontSize',11)
text(0.2,0.2,sprintf('\\beta = %.2f',beta),'FontSize',12,'FontWeight','bold')

xlabel('U_1 (Material uncertainty)')
ylabel('U_2 (Load uncertainty)')

title('Reliability Geometry in Standard Normal Space')

axis equal
axis([-6 6 -6 6])

grid on
set(gca,'FontSize',12,'LineWidth',1.2)

legend('Safe/Failure regions','Failure boundary','Origin','Design point','Reliability direction','Importance samples')

end