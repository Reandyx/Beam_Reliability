function convergence_study(mu,sigma,F,L,I,delta_allow)

% Convergence study for Monte Carlo failure probability
% Evaluates Pf as function of sample size N

% Sample sizes (log-spaced)
N_values = round(logspace(3,6,10));

Pf_values = zeros(length(N_values),1);

for k = 1:length(N_values)

    N = N_values(k);

    % Generate samples
    E_samples = generate_E_samples(mu,sigma,N);

    % Compute deflections
    delta_samples = monte_carlo_deflection(F,L,E_samples,I);

    % Compute failure probability
    [Pf,~] = failure_probability_mc(delta_samples,delta_allow);

    Pf_values(k) = Pf;

end

disp(table(N_values', Pf_values, 'VariableNames', {'N','Pf'}));

% Plot convergence
figure
semilogx(N_values, Pf_values, '-o','LineWidth',1.5)
xlabel('Number of samples N')
ylabel('Failure Probability P_f')
title('Monte Carlo Convergence Study')
grid on

saveas(gcf,'plots/convergence_pf.png')

end