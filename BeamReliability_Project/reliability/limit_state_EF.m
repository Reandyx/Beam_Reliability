function g = limit_state_EF(Z, params)
% LIMIT_STATE_EF
% 2D Limit state in standard normal space
%
% Z      : [2x1] vector [Z_E; Z_F]
% params : struct with fields
%          mu_E, sigma_E
%          mu_F, sigma_F
%          L, I, delta_allow
%
% g(Z) = delta_allow - delta(E,F)
% Failure if g < 0

    % --- Extract Z components ---
    Z_E = Z(1);
    Z_F = Z(2);

    % --- Transform to physical space ---
    E = params.mu_E + params.sigma_E * Z_E;
    F = params.mu_F + params.sigma_F * Z_F;

    % --- Compute deflection via mechanics layer ---
    delta = beam_deflection(F, params.L, E, params.I);

    % --- Limit state ---
    g = params.delta_allow - delta;
end