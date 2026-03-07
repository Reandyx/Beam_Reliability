function plot_sobol_results(S_E,S_F,S_E_total,S_F_total)

values = [
    S_E S_E_total
    S_F S_F_total
];

figure

bar(values)

set(gca,'XTickLabel',{'E','F'})

legend('First Order','Total Order','Location','northwest')

ylabel('Sobol Index')

title('Global Sensitivity Analysis (Sobol Indices)')

grid on

end