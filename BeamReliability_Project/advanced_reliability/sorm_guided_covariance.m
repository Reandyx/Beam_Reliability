function Sigma = sorm_guided_covariance(kappa, dim)
% SORM_GUIDED_COVARIANCE
% Builds covariance matrix for importance sampling using curvature

Sigma = eye(dim);

scale = 1/(1 + abs(kappa));

Sigma = Sigma * scale;

end