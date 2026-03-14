function g = limit_state_E(E, params)

% Limit state function
% g(E) = delta_allow - F*L^3/(3*E*I)

g = params.delta_allow - (params.F .* params.L.^3) ./ (3 .* E .* params.I);

end