function plot_failure_density(params)

N = 50000;

[E,F] = generate_EF_samples(params,N);

g = limit_state_EF(E,F,params);

figure
scatter(E(g>0),F(g>0),10,'b','filled')
hold on
scatter(E(g<0),F(g<0),10,'r','filled')

xlabel('Young Modulus E')
ylabel('Load F')

legend('Safe','Failure')

title('Failure Domain in Parameter Space')

grid on

end