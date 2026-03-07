clear; clc; close all;

%% ===============================
%  TASK 0.4 — Deterministic Baseline
% ================================

F = 1000;        % [N]
L = 1.0;         % [m]
b = 0.05;        % [m]
h = 0.05;        % [m]
E = 210e9;       % [Pa]

I = section_inertia_rect(b,h);
delta = beam_deflection(F,L,E,I);

fprintf('Moment of inertia: %.6e m^4\n', I);
fprintf('Deflection: %.6e m\n', delta);
fprintf('Deflection: %.6f mm\n', delta*1000);

% ---- Small deflection validation ----
delta_ratio = delta / L;
fprintf('delta/L ratio: %.6f\n', delta_ratio);

% Maximum bending stress
M_max = F * L;          % [N*m]
c = h / 2;              % [m]
sigma_max = M_max * c / I;   % [Pa]

sigma_y = 250e6;        % 250 MPa yield steel
safety_factor = sigma_y / sigma_max;

fprintf('Max bending stress: %.2f MPa\n', sigma_max/1e6);
fprintf('Safety factor (yield): %.2f\n', safety_factor);

E_low = 10e9;
delta_lowE = beam_deflection(F,L,E_low,I);

fprintf('Deflection at E = 10 GPa: %.6f m\n', delta_lowE);

%% ===============================
%  TASK 0.6 — Scaling Studies
% ================================

if ~exist('plots','dir')
    mkdir('plots')
end

% ----- Study A: delta vs E -----
E_range = linspace(50e9,250e9,200);
delta_E = beam_deflection(F,L,E_range,I);

figure;
plot(E_range/1e9, delta_E*1000,'LineWidth',2);
xlabel('Young''s Modulus [GPa]');
ylabel('Deflection [mm]');
title('Deflection vs Young''s Modulus');
grid on;
saveas(gcf,'plots/delta_vs_E.png');

% ----- Study B: delta vs L -----
L_range = linspace(0.5,2.0,200);
delta_L = beam_deflection(F,L_range,E,I);

figure;
plot(L_range, delta_L*1000,'LineWidth',2);
xlabel('Length [m]');
ylabel('Deflection [mm]');
title('Deflection vs Length');
grid on;
saveas(gcf,'plots/delta_vs_L.png');

% ----- Study C: delta vs F -----
F_range = linspace(100,2000,200);
delta_F = beam_deflection(F_range,L,E,I);

figure;
plot(F_range, delta_F*1000,'LineWidth',2);
xlabel('Force [N]');
ylabel('Deflection [mm]');
title('Deflection vs Force');
grid on;
saveas(gcf,'plots/delta_vs_F.png');

%% ===============================
%  TASK 0.7 — Sensitivity Verification
% ================================

h_fd = 1e6;  % finite difference step [Pa]

analytical = analytical_sensitivity(F,L,E,I);

delta_plus  = beam_deflection(F,L,E+h_fd,I);
delta_minus = beam_deflection(F,L,E-h_fd,I);

numerical = (delta_plus - delta_minus) / (2*h_fd);

error_rel = abs(analytical - numerical) / abs(analytical);

fprintf('\nSensitivity verification:\n');
fprintf('Analytical: %.6e\n', analytical);
fprintf('Numerical : %.6e\n', numerical);
fprintf('Relative error: %.6e\n', error_rel);

%% ===============================
%  TASK 0.8 — Non-Dimensional Check
% ================================

delta_bar = (3*E*I*delta) / (F*L^3);

fprintf('\nNon-dimensional consistency check:\n');
fprintf('|delta_bar - 1| = %.6e\n', abs(delta_bar - 1));