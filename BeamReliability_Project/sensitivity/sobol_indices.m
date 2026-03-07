function [S_E,S_F,S_E_total,S_F_total] = sobol_indices(params,N)

% Number of variables
d = 2;

%% Generate Sobol base matrices

A = randn(N,d);
B = randn(N,d);

% Transform to physical variables
E_A = params.mu_E + params.sigma_E*A(:,1);
F_A = params.mu_F + params.sigma_F*A(:,2);

E_B = params.mu_E + params.sigma_E*B(:,1);
F_B = params.mu_F + params.sigma_F*B(:,2);

% Mixed matrices
E_AB1 = E_B;
F_AB1 = F_A;

E_AB2 = E_A;
F_AB2 = F_B;

%% Model evaluation

YA = beam_deflection(F_A,params.L,E_A,params.I);
YB = beam_deflection(F_B,params.L,E_B,params.I);

YAB1 = beam_deflection(F_AB1,params.L,E_AB1,params.I);
YAB2 = beam_deflection(F_AB2,params.L,E_AB2,params.I);

%% Variance

VarY = var([YA;YB]);

%% First order indices

S_E = mean(YB .* (YAB1 - YA)) / VarY;
S_F = mean(YB .* (YAB2 - YA)) / VarY;

%% Total order indices

S_E_total = mean((YA - YAB1).^2) / (2*VarY);
S_F_total = mean((YA - YAB2).^2) / (2*VarY);

end