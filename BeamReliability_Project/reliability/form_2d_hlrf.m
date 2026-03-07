function [beta, Z_star, Pf_FORM, Z_history] = form_2d_hlrf(params)

    max_iter = 100;
    tol = 1e-6;

    Z = [0; 0];
    Z_history = Z';
    for k = 1:max_iter

        g = limit_state_EF(Z, params);
        grad = gradient_EF(Z, params);

        norm_grad = norm(grad);

        if norm_grad < 1e-12
            error('Gradient too small in HL-RF.');
        end

        % Standard HL-RF update
        Z_new = (grad' * Z - g) / (norm_grad^2) * grad;
        
        Z_history = [Z_history; Z_new'];

        if isfield(params,'verbose') && params.verbose
            fprintf('Iter %d | ||Z|| = %.6f\n', k, norm(Z_new));
        end

        if norm(Z_new - Z) < tol
            Z = Z_new;
            break;
        end

        Z = Z_new;

    end

    Z_star = Z;
    beta = norm(Z_star);
    Pf_FORM = normcdf(-beta);
    
    if k == max_iter
    warning("HL-RF reached maximum iterations without convergence.");
    end
    
end