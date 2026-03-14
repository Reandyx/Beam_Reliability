function ESS = compute_ess(weights)
% COMPUTE_ESS
% Computes effective sample size of importance sampling weights
%
% ESS = (sum(w)^2) / sum(w^2)

w = weights(:);

ESS = (sum(w)^2) / sum(w.^2);

end