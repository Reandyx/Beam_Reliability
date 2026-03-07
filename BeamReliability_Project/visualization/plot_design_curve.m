function plot_design_curve(h_vals, beta_vals, h_opt, beta_opt, beta_target)

figure
plot(h_vals, beta_vals,'LineWidth',2)
hold on

yline(beta_target,'r--','LineWidth',2)

plot(h_opt,beta_opt,'ko','MarkerSize',10,'MarkerFaceColor','k')

xlabel('Beam Height h [m]')
ylabel('Reliability Index \beta')

title('Reliability-Based Design Curve')

legend('\beta(h)','Target Reliability','Optimal Design','Location','NorthWest')

grid on
box on

end