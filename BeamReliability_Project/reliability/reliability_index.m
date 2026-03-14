function beta_analytical = reliability_index(params)

E_star = (params.F * params.L^3) / ...
         (3 * params.I * params.delta_allow);

Z_star = (E_star - params.mu_E) / params.sigma_E;

beta_analytical = abs(Z_star);

end