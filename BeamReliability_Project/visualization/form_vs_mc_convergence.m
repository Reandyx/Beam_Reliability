function form_vs_mc_convergence(params)

rng(1)

sample_sizes = round(logspace(3,6,8));

Pf_MC = zeros(size(sample_sizes));
CI_MC = zeros(size(sample_sizes));

[~,~,Pf_FORM] = form_2d_hlrf(params);

for i = 1:length(sample_sizes)

    N = sample_sizes(i);

    [E,F] = generate_EF_samples( ...
        params.mu_E,params.sigma_E, ...
        params.mu_F,params.sigma_F,N);

    [Pf_MC(i),CI_MC(i)] = failure_probability_mc_2d(E,F,params);

    fprintf("N = %d, Pf_MC = %.6e\n",N,Pf_MC(i))

end


%% ------------------------------------------------------------
% Prepare data for log-scale plotting
%% ------------------------------------------------------------

eps_val = 1e-12;

Pf_plot = Pf_MC;
Pf_plot(Pf_plot==0) = eps_val;

lower = max(Pf_MC - CI_MC,0);
upper = Pf_MC + CI_MC;

lower(lower==0) = eps_val;
upper(upper==0) = eps_val;


%% ------------------------------------------------------------
% Plot
%% ------------------------------------------------------------

figure
hold on
set(gcf,'Color','w')

fill([sample_sizes fliplr(sample_sizes)], ...
     [lower fliplr(upper)], ...
     [0.75 0.85 1], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.25)

semilogx(sample_sizes,Pf_plot,'o-', ...
         'LineWidth',1.8, ...
         'MarkerSize',7, ...
         'Color',[0 0.45 0.74])

yline(Pf_FORM,'r--','LineWidth',3)

set(gca,'YScale','log')
set(gca,'XScale','log')

ylim([1e-16 1e-5])

xlabel('Number of MC samples N')
ylabel('Failure probability P_f')

title('FORM vs Monte Carlo convergence')

text(2e3,1e-8,...
'No failures observed (rare-event regime)',...
'FontSize',11)

legend('95% MC confidence interval', ...
       'Monte Carlo estimate', ...
       'FORM prediction', ...
       'Location','best')

grid on
box on
set(gca,'FontSize',12,'LineWidth',1)

end