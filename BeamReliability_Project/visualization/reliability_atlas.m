function reliability_atlas(params)

% reasonable resolution
cov_E = linspace(0.01,0.30,100);
cov_F = linspace(0.01,0.30,100);

[CE,CF] = meshgrid(cov_E,cov_F);

Pf_map = NaN(size(CE));

fprintf('Computing reliability atlas...\n')

for i = 1:size(CE,1)
for j = 1:size(CE,2)

    params.sigma_E = CE(i,j) * params.mu_E;
    params.sigma_F = CF(i,j) * params.mu_F;

    try
        [~,~,Pf] = form_2d_hlrf(params);

        if ~isnan(Pf) && Pf > 0
            Pf_map(i,j) = Pf;
        end

    catch
        Pf_map(i,j) = NaN;
    end

end
end

% fill missing values AFTER the loop
Pf_map = fillmissing(Pf_map,'nearest');

% smooth the field
Pf_map = imgaussfilt(Pf_map,1.2);

% clamp extreme values
Pf_map(Pf_map < 1e-12) = 1e-12;

logPf = log10(Pf_map);

figure
hold on

imagesc(cov_E,cov_F,logPf)
set(gca,'YDir','normal')

colormap(turbo)
colorbar

% readable probability range
caxis([-10 -2])

% reliability contours
contour(CE,CF,logPf,[-10 -8 -6 -4 -3],'k','LineWidth',1)

xlabel('CoV_E (Material Uncertainty)')
ylabel('CoV_F (Load Uncertainty)')
title('Reliability Atlas: log_{10}(P_f)')

axis square
grid on

set(gca,'FontSize',12,'LineWidth',1.2)

end