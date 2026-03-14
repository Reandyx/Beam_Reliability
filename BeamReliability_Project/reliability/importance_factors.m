function alpha = importance_factors(Z_star, beta)
% IMPORTANCE_FACTORS
% Computes FORM importance factors

    if beta == 0
        error('Beta is zero, cannot compute importance factors.');
    end

    alpha = Z_star / beta;
    alpha = alpha / norm(alpha);
    
end