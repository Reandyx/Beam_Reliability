function plot_sobol_results(S_E,S_F,S_E_total,S_F_total)

data = [
    S_E      S_E_total
    S_F      S_F_total
];

figure
hold on

b = bar(data,'grouped');

b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];

xticks([1 2])
xticklabels({'E','F'})
set(gca,'XTickLabel',{'E','F'})

ylabel('Sobol Index')
xlabel('Input variable')

title('Global Sensitivity Analysis of Beam Reliability (Sobol Indices)')

legend({'First Order','Total Order'},'Location','northwest')

ylim([0 max(data(:))*1.15])

grid on
set(gca,'YGrid','on','XGrid','off')

set(gca,'FontSize',12)

%% Add labels

for i = 1:2
    xtips = b(i).XEndPoints;
    ytips = b(i).YEndPoints;
    labels = compose('%.3f',b(i).YData);
    text(xtips,ytips,labels,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontWeight','bold');
end

end