function Pf_SORM = sorm_correction(params)

% Compute FORM design point
[beta, Z_star, Pf_FORM] = form_2d_hlrf(params);

grad = gradient_EF(Z_star, params);
grad_norm = norm(grad);

% Numerical Hessian
h = 1e-5;
H = zeros(2);

for i = 1:2
    for j = 1:2

        Z1 = Z_star;
        Z2 = Z_star;
        Z3 = Z_star;
        Z4 = Z_star;

        Z1(i)=Z1(i)+h; Z1(j)=Z1(j)+h;
        Z2(i)=Z2(i)+h; Z2(j)=Z2(j)-h;
        Z3(i)=Z3(i)-h; Z3(j)=Z3(j)+h;
        Z4(i)=Z4(i)-h; Z4(j)=Z4(j)-h;

        g1 = limit_state_EF(Z1,params);
        g2 = limit_state_EF(Z2,params);
        g3 = limit_state_EF(Z3,params);
        g4 = limit_state_EF(Z4,params);

        H(i,j) = (g1 - g2 - g3 + g4) / (4*h^2);

    end
end

% Curvature approximation
kappa = (grad' * H * grad) / (grad_norm^3);

% SORM correction
Pf_SORM = normcdf(-beta) / sqrt(1 + beta * kappa);

end