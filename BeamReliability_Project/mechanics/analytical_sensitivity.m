function ddelta_dE = analytical_sensitivity(F,L,E,I)
%ANALYTICAL_SENSITIVITY Analytical derivative d(delta)/dE
%
%   d(delta)/dE = -F*L^3 / (3*I*E^2)

ddelta_dE = -(F .* L.^3) ./ (3 .* I .* E.^2);

end