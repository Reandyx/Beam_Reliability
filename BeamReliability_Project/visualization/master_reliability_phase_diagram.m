function master_reliability_phase_diagram(~)

cov_E = linspace(0,0.30,300);
cov_F = linspace(0,0.30,300);

[CE,CF] = meshgrid(cov_E,cov_F);

beta = 1 ./ sqrt(CE.^2 + CF.^2 + 1e-6);

figure('Color','w')
hold on

% --- REGION SHADING (high contrast)
patch([0 0.3 0.3 0],[0 0 0.15 0.15],[0.88 0.93 1],'EdgeColor','none')
patch([0 0.3 0.3 0],[0.15 0.15 0.3 0.3],[1 0.88 0.88],'EdgeColor','none')

% --- RELIABILITY CONTOURS
[C,h] = contour(CE,CF,beta,[3 4 5 6],'k','LineWidth',2);
clabel(C,h,'FontSize',11,'Color','k')

% --- FORM VALIDITY BOUNDARY
plot([0 0.30],[0.15 0.15],'k','LineWidth',3)

% --- DOMINANCE LINE
plot([0 0.30],[0 0.30],'k--','LineWidth',2)

% --- LABELS
text(0.03,0.035,'FORM accurate','FontSize',12,'FontWeight','bold')

text(0.04,0.23,'Nonlinear regime (SORM)','FontSize',12,'FontWeight','bold')

text(0.015,0.18,'Load dominated','Rotation',45,'FontSize',11)

text(0.22,0.05,'Material dominated','Rotation',45,'FontSize',11)

% --- AXES
xlabel('CoV_E  (Material uncertainty)','FontSize',12)
ylabel('CoV_F  (Load uncertainty)','FontSize',12)

title('Reliability Phase Diagram','FontSize',14)

axis([0 0.30 0 0.30])
axis square
grid on

set(gca,'FontSize',12,'LineWidth',1.3)

end