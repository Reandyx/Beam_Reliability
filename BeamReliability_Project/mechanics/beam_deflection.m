function delta = beam_deflection(F,L,E,I)
%BEAM_DEFLECTION Computes tip deflection of a cantilever beam
%   delta = F*L^3 / (3*E*I)
%   Inputs:
%       F - Force [N]
%       L - Length [m]
%       E - Young's modulus [Pa = N/m^2]
%       I - Second moment of area [m^4]
%   Output:
%       delta - Deflection [m]
%   Dimensional Verification:
%       F   -> [N]
%       L^3 -> [m^3]
%       E   -> [N/m^2]
%       I   -> [m^4]
%       Numerator: N * m^3
%       Denominator: (N/m^2) * m^4 = N * m^2
%       Final unit:
%           (N*m^3)/(N*m^2) = m
%       Result in meters

delta = (F .* L.^3) ./ (3 .* E .* I);

end