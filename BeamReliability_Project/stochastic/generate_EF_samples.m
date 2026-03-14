function [E_samples, F_samples] = generate_EF_samples(mu_E, sigma_E, ...
                                                      mu_F, sigma_F, N)

    E_samples = mu_E + sigma_E .* randn(N,1);
    F_samples = mu_F + sigma_F .* randn(N,1);
    
    % Remove nonphysical values
    E_samples(E_samples <= 0) = mu_E;
    F_samples(F_samples <= 0) = mu_F;
end