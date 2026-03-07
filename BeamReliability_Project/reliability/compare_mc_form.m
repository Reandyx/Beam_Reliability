function compare_mc_form(Pf_MC, mu_E, sigma_E, F, L, I, delta_allow)

[beta_solver, Pf_FORM, Z_star, E_star] = ...
    form_1d_solver(mu_E, sigma_E, F, L, I, delta_allow);

beta_analytical = reliability_index(mu_E, sigma_E, F, L, I, delta_allow);

abs_beta_error = abs(beta_solver - beta_analytical);
rel_beta_error = abs_beta_error / beta_analytical;

abs_Pf_error = abs(Pf_MC - Pf_FORM);
rel_Pf_error = abs_Pf_error / Pf_MC;

fprintf('\n=== FORM vs Monte Carlo Comparison ===\n');
fprintf('Design point E*       = %.6f\n', E_star);
fprintf('Reliability index β   = %.6f\n', beta_solver);
fprintf('Pf_FORM               = %.6e\n', Pf_FORM);
fprintf('Pf_MC                 = %.6e\n', Pf_MC);
fprintf('Absolute Pf error     = %.6e\n', abs_Pf_error);
fprintf('Relative Pf error     = %.6e\n', rel_Pf_error);
fprintf('β analytical          = %.6f\n', beta_analytical);
fprintf('β solver difference   = %.6e\n', abs_beta_error);
fprintf('β relative difference = %.6e\n', rel_beta_error);

end