function delta_samples = monte_carlo_deflection(F,L,E_samples,I)

% Computes deflection samples using deterministic beam model
% No reimplementation of mechanics allowed
% Output: column vector

delta_samples = beam_deflection(F,L,E_samples,I);

delta_samples = delta_samples(:); % enforce column format

end