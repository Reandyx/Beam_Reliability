function grad_g = gradient_EF(Z, params)
% GRADIENT_EF
% Analytical gradient of limit state in Z-space
%
% Returns:
% grad_g : [2x1] vector ∇g(Z)

    % --- Extract Z ---
    Z_E = Z(1);
    Z_F = Z(2);

    % --- Transform to physical space ---
    E = params.mu_E + params.sigma_E * Z_E;
    F = params.mu_F + params.sigma_F * Z_F;

    L = params.L;
    I = params.I;

    % --- Physical derivatives ---
    dgdE = (F * L^3) / (3 * I * E^2);
    dgdF = - (L^3) / (3 * E * I);

    % --- Chain rule to Z-space ---
    dgdZ_E = dgdE * params.sigma_E;
    dgdZ_F = dgdF * params.sigma_F;

    grad_g = [dgdZ_E; dgdZ_F];

end