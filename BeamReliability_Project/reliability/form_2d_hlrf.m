function [beta, Z_star, Pf_FORM, Z_history] = form_2d_hlrf(params)

max_iter = 100;
tol = 1e-6;

Z = [0;0];
Z_history = Z';

for k = 1:max_iter

    g = limit_state_EF(Z,params);
    grad = gradient_EF(Z,params);

    norm_grad = norm(grad);

    if norm_grad < 1e-12
        beta = NaN;
        Z_star = Z;
        Pf_FORM = NaN;
        return
    end

    Z_proj = ((grad' * Z) - g) / (norm_grad^2) * grad;

    lambda = 0.5;          % damping factor
    Z_new = Z + lambda*(Z_proj - Z);
    
    Z_history = [Z_history; Z_new'];

    if norm(Z_new - Z) < tol
        Z = Z_new;
        break
    end

    Z = Z_new;

end

Z_star = Z;

beta = norm(Z_star);

Pf_FORM = normcdf(-beta);
disp(beta)
disp(Z_star)
end