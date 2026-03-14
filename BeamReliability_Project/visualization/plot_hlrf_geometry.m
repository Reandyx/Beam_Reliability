function plot_hlrf_geometry(U_history, design_point, params)

u = -4:0.1:4;
[UE,UF] = meshgrid(u,u);

g = zeros(size(UE));

for i=1:numel(UE)

    Z = [UE(i);UF(i)];
    g(i) = limit_state_EF(Z,params);

end

figure
hold on
grid on
box on

% --- Failure surface ---
[C,h] = contour(UE,UF,g,[0 0],'LineWidth',3,'Color',[0 0.8 0.8]);

% --- Failure region shading ---
contourf(UE,UF,g,[-100 0],'FaceColor',[1 0.3 0.3],'FaceAlpha',0.15,'LineStyle','none')

% --- Safe region shading ---
contourf(UE,UF,g,[0 100],'FaceColor',[0.3 0.7 1],'FaceAlpha',0.15,'LineStyle','none')

% --- HL-RF trajectory ---
plot(U_history(:,1),U_history(:,2),'o-','LineWidth',2,'Color',[0.2 0.6 1],'MarkerSize',7)

% --- Design point ---
plot(design_point(1),design_point(2),'ro','MarkerSize',12,'MarkerFaceColor','r')

% --- Origin ---
plot(0,0,'ko','MarkerFaceColor','k')

% --- Reliability index vector ---
quiver(0,0,design_point(1),design_point(2),0,...
'LineWidth',2,'Color','k','MaxHeadSize',0.5)

% --- Beta label ---
beta = norm(design_point);
text(design_point(1)/2, design_point(2)/2,...
    sprintf('\\beta = %.2f',beta),...
    'FontSize',12,'FontWeight','bold')

% --- Tangent line (FORM approximation) ---
grad = gradient_EF(design_point,params);

t = -4:0.1:4;

tangent = -(grad(1)/grad(2))*(t-design_point(1)) + design_point(2);

plot(t,tangent,'--','LineWidth',2,'Color',[1 0.6 0])

% --- Normal vector to failure surface (gradient) ---

grad = gradient_EF(design_point,params);

grad_unit = grad / norm(grad);

scale = 0.8;

quiver(design_point(1),design_point(2),...
       grad_unit(1)*scale,...
       grad_unit(2)*scale,...
       0,...
       'Color',[0 1 0],...
       'LineWidth',2,...
       'MaxHeadSize',0.5)

xlabel('U_E')
ylabel('U_F')

title('FORM Reliability Geometry')

legend({'Failure Surface g(U)=0',...
        'Failure Region',...
        'Safe Region',...
        'HL-RF Trajectory',...
        'Design Point U*',...
        'Origin',...
        'Reliability Vector',...
        'FORM Tangent Approximation',...
        'Failure Surface Normal ∇g'},...
        'Location','northwest')

axis equal

end