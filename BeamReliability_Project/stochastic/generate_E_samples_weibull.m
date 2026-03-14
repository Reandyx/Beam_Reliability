function E = generate_E_samples_weibull(mu,sigma,N)

k = (sigma/mu)^(-1.086);

lambda = mu / gamma(1 + 1/k);

E = wblrnd(lambda,k,N,1);

end