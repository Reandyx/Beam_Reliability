function E = generate_E_samples_lognormal(mu,sigma,N)

phi = sqrt(log(1 + (sigma/mu)^2));
mu_ln = log(mu) - 0.5*phi^2;

E = lognrnd(mu_ln,phi,N,1);

end