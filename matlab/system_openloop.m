%% CPEN 401 -- Group 16: Liquid Flow Control in a Pipeline Pumping Station
%% Script 01: System Definition & Open-Loop Analysis
%  Run this script first. It defines the plant and generates:
%    fig01_openloop_step.png
%    fig02_pole_zero.png
%    fig03_bode_openloop.png
%    fig04_root_locus.png
%  Save figures to the 'figures/' subfolder.

clear; clc; close all;

%% ── SYSTEM PARAMETERS ──────────────────────────────────────────────────────
rho   = 1000;       % fluid density [kg/m^3]
mu    = 0.001;      % dynamic viscosity [Pa.s]
L     = 100;        % pipe length [m]
D     = 0.10;       % pipe diameter [m]
A_c   = pi*D^2/4;   % cross-section area [m^2]
Ih    = rho*L/A_c;  % hydraulic inductance [kg/m^4]
Rp    = 128*mu*L/(pi*D^4);  % Hagen-Poiseuille resistance [Pa.s/m^3]
Rv    = 500;        % valve resistance [Pa.s/m^3]
Rtot  = Rp + Rv;
tau_p = 2.0;        % pump time constant [s]
Kpump = 50000;      % pump gain [Pa/unit]

% Normalised model poles
T1 = 2.0;    % pump time constant [s]
T2 = 7.05;   % pipeline hydraulic time constant [s]
p1 = -1/T1;  % -0.5000 rad/s
p2 = -1/T2;  % -0.1418 rad/s

fprintf('=== SYSTEM PARAMETERS ===\n');
fprintf('Pipe area       A  = %.4e m^2\n', A_c);
fprintf('Hydr. inductance Ih = %.2f kg/m^4\n', Ih);
fprintf('Pipe resistance  Rp = %.1f Pa.s/m^3\n', Rp);
fprintf('Total resistance Rt = %.1f Pa.s/m^3\n', Rtot);
fprintf('Pole p1 = %.4f rad/s  (pump)\n', p1);
fprintf('Pole p2 = %.4f rad/s  (pipeline)\n', p2);

%% ── PLANT TRANSFER FUNCTION ────────────────────────────────────────────────
% G(s) = K / [(T1*s+1)(T2*s+1)]   K=1 (normalised)
num_G = [abs(p1)*abs(p2)];                    % 0.07092
den_G = [1, -(p1+p2), p1*p2];                % [1, 0.6418, 0.07092]
G = tf(num_G, den_G);
fprintf('\nPlant G(s):\n'); display(G);

%% ── FIG 01: OPEN-LOOP STEP RESPONSE ───────────────────────────────────────
figure('Name','Open-Loop Step','Position',[100 100 800 420]);
[y,t] = step(G, 0:0.1:200);
plot(t, y, 'b', 'LineWidth', 2); hold on;
yss = dcgain(G);
plot([0 200], [yss yss], 'r--', 'LineWidth', 1.2);
xlabel('Time (s)','FontSize',12);
ylabel('Flow Rate (normalised)','FontSize',12);
title('Open-Loop Step Response -- Pipeline Pumping Station','FontSize',13);
legend('y(t)', sprintf('Steady-state = %.3f', yss), 'Location','southeast');
grid on; box on;
saveas(gcf, 'figures/fig01_openloop_step.png');
fprintf('Saved fig01_openloop_step.png\n');

%% ── FIG 02: POLE-ZERO MAP ─────────────────────────────────────────────────
figure('Name','Pole-Zero Map','Position',[100 100 600 500]);
pzmap(G);
title('Pole-Zero Map -- Open-Loop Plant G(s)','FontSize',13);
grid on;
% Annotate poles
hold on;
poles_G = roots(den_G);
for k = 1:length(poles_G)
    text(real(poles_G(k))+0.005, imag(poles_G(k))+0.01, ...
         sprintf('  p = %.4f', real(poles_G(k))), 'FontSize', 10);
end
saveas(gcf, 'figures/fig02_pole_zero.png');
fprintf('Saved fig02_pole_zero.png\n');

%% ── FIG 03: BODE PLOT (OPEN LOOP) ─────────────────────────────────────────
figure('Name','Bode Plot','Position',[100 100 820 600]);
bode(G);
title('Bode Plot -- Open-Loop Plant G(s)','FontSize',13);
grid on;
% Add frequency markers at pole locations
[gm, pm, wpc, wgc] = margin(G);
fprintf('\nOpen-loop gain margin  = %.2f dB\n', 20*log10(gm));
fprintf('Open-loop phase margin = %.2f deg\n', pm);
saveas(gcf, 'figures/fig03_bode_openloop.png');
fprintf('Saved fig03_bode_openloop.png\n');

%% ── FIG 04: ROOT LOCUS ────────────────────────────────────────────────────
figure('Name','Root Locus','Position',[100 100 700 600]);
rlocus(G);
title('Root Locus -- Plant G(s) with Unity Feedback','FontSize',13);
sgrid;
xlim([-3 0.5]); ylim([-1.5 1.5]);
saveas(gcf, 'figures/fig04_root_locus.png');
fprintf('Saved fig04_root_locus.png\n');

fprintf('\nScript 01 complete. All open-loop figures saved.\n');
