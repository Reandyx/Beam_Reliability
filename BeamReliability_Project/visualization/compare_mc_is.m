function compare_mc_is(mc_estimates, is_estimates, sample_sizes)

figure
semilogx(sample_sizes, mc_estimates, 'r-o', 'LineWidth',2)
hold on
semilogx(sample_sizes, is_estimates, 'b-s', 'LineWidth',2)

xlabel('Number of Samples')
ylabel('Estimated Failure Probability')
title('Monte Carlo vs Importance Sampling Convergence')

legend('Monte Carlo','Importance Sampling')
grid on

end