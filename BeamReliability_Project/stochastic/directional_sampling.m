function U = directional_sampling(U_star, N)
% DIRECTIONAL_SAMPLING
% Generates samples aligned with FORM reliability direction

dim = length(U_star);

alpha = U_star / norm(U_star);

% radial component along reliability direction
r = randn(N,1);

% orthogonal noise
V = randn(N,dim);

% remove projection along alpha
V = V - (V*alpha)*alpha';

U = r*alpha' + V;

end