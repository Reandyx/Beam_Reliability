function E_samples = generate_E_samples(mu,sigma,N)

% Generates N Gaussian samples of Young's modulus
% E ~ N(mu, sigma)
% Output: column vector [N x 1]
% Units: Pa (SI)

E_samples = mu + sigma .* randn(N,1);

% Remove nonphysical values
E_samples(E_samples <= 0) = mu;

end