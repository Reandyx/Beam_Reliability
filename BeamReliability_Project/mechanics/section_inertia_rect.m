function I = section_inertia_rect(b,h)
%SECTION_INERTIA_RECT Computes second moment of area for rectangle
%   I = (b*h^3)/12
%   Inputs:
%       b - width  [m]
%       h - height [m]
%   Output:
%       I - second moment of inertia [m^4]
%   Unit verification:
%       b  -> [m]
%       h^3 -> [m^3]
%       I = m * m^3 = m^4

I = (b .* h.^3) ./ 12;

end