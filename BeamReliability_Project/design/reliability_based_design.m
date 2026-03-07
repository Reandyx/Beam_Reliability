function [h_opt, beta_opt, I_opt] = reliability_based_design(params, beta_target, b)

% RELIABILITY_BASED_DESIGN
% Robust RBD solver using FORM + root finding

% search interval
h_min = 0.01;
h_max = 0.30;

f = @(h) beta_difference(h, params, beta_target, b);

% Check root bracket
f_min = f(h_min);
f_max = f(h_max);

if sign(f_min) == sign(f_max)
    error('Root not bracketed. Increase search interval.');
end

% Solve for optimal height
h_opt = fzero(f,[h_min h_max]);

% Final inertia
I_opt = b*h_opt^3/12;

params.I = I_opt;

[beta_opt,~,~] = form_2d_hlrf(params);

end


function val = beta_difference(h, params, beta_target, b)

I = b*h^3/12;
params.I = I;

try
    [beta,~,~] = form_2d_hlrf(params);
catch
    beta = -10;
end

val = beta - beta_target;

end