function [U_history, design_point, curvature] = hlrf_trajectory_analysis(params)

% Run FORM solver
[beta, Z_star, Pf_FORM, Z_history] = form_2d_hlrf(params);

U_history = Z_history;
design_point = Z_star;

% --- Curvature estimation via Hessian ---

Z = Z_star;

h = 1e-4;

% second derivatives (finite difference)

g0 = limit_state_EF(Z, params);

Z1 = Z + [h;0];
Z2 = Z - [h;0];
g_EE = (limit_state_EF(Z1,params) - 2*g0 + limit_state_EF(Z2,params)) / h^2;

Z1 = Z + [0;h];
Z2 = Z - [0;h];
g_FF = (limit_state_EF(Z1,params) - 2*g0 + limit_state_EF(Z2,params)) / h^2;

Z1 = Z + [h;h];
Z2 = Z + [h;-h];
Z3 = Z + [-h;h];
Z4 = Z + [-h;-h];

g_EF = (limit_state_EF(Z1,params) - limit_state_EF(Z2,params) ...
       -limit_state_EF(Z3,params) + limit_state_EF(Z4,params)) / (4*h^2);

H = [g_EE g_EF;
     g_EF g_FF];

curvature = norm(H);

end