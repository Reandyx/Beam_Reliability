function [h_opt,beta_opt,I_opt] = reliability_based_design(params,beta_target,b)

    beta_fun = @(h) beta_at_h(h,params,b) - beta_target;

    % search interval
    h_low  = 0.08;
    h_high = 0.25;

    % solve β(h)=β_target
    h_opt = fzero(beta_fun,[h_low h_high]);

    % compute inertia
    I_opt = b*h_opt^3/12;

    params.h = h_opt;
    params.I = I_opt;

    [beta_opt,~,~] = form_2d_hlrf(params);

end

function beta = beta_at_h(h,params,b)

    params.h = h;
    params.I = b*h^3/12;

    try
        [beta,~,~] = form_2d_hlrf(params);

        if ~isfinite(beta)
            beta = -10;   % treat failure as very unsafe
        end

    catch
        beta = -10;       % also treat solver crash as unsafe
    end

end