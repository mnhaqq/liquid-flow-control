%% CPEN 401 -- Group 16: Liquid Flow Control in a Pipeline Pumping Station
%% Script 02: Classical Controller Design (P, PI, PD, PID, Lead, Lag, Lead-Lag)
%  Generates:
%    fig05_step_comparison.png
%    fig06_bode_controllers.png
%    fig07_disturbance.png
%    fig08_sensitivity.png
%    fig09_noise.png
%    fig10_saturation.png

clear; clc; close all;

%% ── PLANT ──────────────────────────────────────────────────────────────────
T1 = 2.0; T2 = 7.05;
p1 = -1/T1; p2 = -1/T2;
num_G = [abs(p1)*abs(p2)];
den_G = [1, -(p1+p2), p1*p2];
G = tf(num_G, den_G);
t_sim = 0:0.1:200;

%% ── CONTROLLER DEFINITIONS ─────────────────────────────────────────────────
% P Controller
Kp_P = 5.0;
C_P  = tf(Kp_P, 1);

% PI Controller:  C(s) = Kp*(Ti*s+1)/(Ti*s)
Kp_PI = 6.0; Ti_PI = 20.0;
C_PI  = tf(Kp_PI*[Ti_PI, 1], [Ti_PI, 0]);

% PD Controller:  C(s) = Kp*(Td*s+1)
Kp_PD = 5.0; Td_PD = 1.5;
C_PD  = tf(Kp_PD*[Td_PD, 1], 1);

% PID Controller: C(s) = (Kd*s^2 + Kp*s + Ki)/s
Kp_PID = 8.0; Ki_PID = 0.4; Kd_PID = 3.0;
C_PID  = tf([Kd_PID, Kp_PID, Ki_PID], [1, 0]);

% Lead Compensator: C(s) = K*(s+z)/(s+p)
K_lead = 12.0; z_lead = 0.15; p_lead = 0.75;
C_lead = tf(K_lead*[1, z_lead], [1, p_lead]);

% Lag Compensator: C(s) = K*(s+z)/(s+p),  z > p
K_lag = 6.5; z_lag = 0.08; p_lag = 0.008;
C_lag = tf(K_lag*[1, z_lag], [1, p_lag]);

% Lead-Lag Compensator
K_ll = 10.0;
C_ll = tf(K_ll*conv([1, z_lead],[1, z_lag]), conv([1, p_lead],[1, p_lag]));

% Store all controllers with labels
controllers = {C_P, C_PI, C_PD, C_PID, C_lead, C_lag, C_ll};
ctrl_names  = {'P', 'PI', 'PD', 'PID', 'Lead', 'Lag', 'Lead-Lag'};
colors = {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.2 0.7 0.3], ...
          [0.7 0.1 0.7], [1.0 0.5 0.0], [0.0 0.6 0.6], [0.5 0.3 0.0]};
styles = {'-','-','-','--','-','-','-.'};

%% ── CLOSED-LOOP SYSTEMS ────────────────────────────────────────────────────
CL = cell(1, length(controllers));
for k = 1:length(controllers)
    CL{k} = feedback(controllers{k}*G, 1);
end

%% ── STABILITY MARGINS ──────────────────────────────────────────────────────
fprintf('%-12s  %10s  %10s  %10s  %10s\n', 'Controller','GM (dB)','PM (deg)','Rise T (s)','Settle T (s)');
fprintf('%s\n', repmat('-',1,60));
for k = 1:length(controllers)
    L = controllers{k}*G;
    [gm,pm] = margin(L);
    gm_dB = 20*log10(gm);
    if isinf(gm_dB), gm_str = 'Inf'; else, gm_str = sprintf('%.1f', gm_dB); end
    [y,t] = step(CL{k}, t_sim);
    % Rise time
    yss = y(end);
    idx10 = find(y >= 0.10*yss, 1); idx90 = find(y >= 0.90*yss, 1);
    if ~isempty(idx10) && ~isempty(idx90)
        tr = t(idx90) - t(idx10);
    else, tr = NaN; end
    % Settling time
    band = 0.02*yss;
    idx_out = find(abs(y - yss) > band);
    if ~isempty(idx_out), ts = t(idx_out(end)); else, ts = 0; end
    fprintf('%-12s  %10s  %10.1f  %10.1f  %10.1f\n', ctrl_names{k}, gm_str, pm, tr, ts);
end

%% ── FIG 05: STEP RESPONSE COMPARISON ──────────────────────────────────────
figure('Name','Step Comparison','Position',[100 100 900 540]);
hold on;
for k = 1:length(controllers)
    [y,t] = step(CL{k}, t_sim);
    plot(t, y, 'Color', colors{k}, 'LineStyle', styles{k}, 'LineWidth', 1.8, ...
         'DisplayName', ctrl_names{k});
end
plot(t_sim, ones(size(t_sim)), 'k:', 'LineWidth', 1.2, 'DisplayName', 'Reference');
hold off;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Normalised Flow Rate', 'FontSize', 12);
title('Closed-Loop Step Response -- Controller Comparison', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 10);
ylim([-0.05 1.5]); xlim([0 200]);
grid on; box on;
saveas(gcf, 'figures/fig05_step_comparison.png');
fprintf('\nSaved fig05_step_comparison.png\n');

%% ── FIG 06: BODE PLOTS (CONTROLLERS) ──────────────────────────────────────
selected = {C_PI, C_PID, C_lead, C_ll};
sel_names = {'PI', 'PID', 'Lead', 'Lead-Lag'};
sel_colors = {[0.8 0.2 0.2],[0.7 0.1 0.7],[1.0 0.5 0.0],[0.5 0.3 0.0]};

figure('Name','Bode Controllers','Position',[100 100 860 640]);
w = logspace(-3, 2, 2000);
for k = 1:length(selected)
    L = selected{k}*G;
    [mag, ph] = bode(L, w);
    mag = squeeze(mag); ph = squeeze(ph);
    subplot(2,1,1); hold on;
    semilogx(w, 20*log10(mag), 'Color', sel_colors{k}, 'LineWidth', 1.8, ...
             'DisplayName', sel_names{k});
    subplot(2,1,2); hold on;
    semilogx(w, ph, 'Color', sel_colors{k}, 'LineWidth', 1.8, ...
             'DisplayName', sel_names{k});
end
subplot(2,1,1);
yline(0, 'k--', 'LineWidth', 1, 'Label', '0 dB');
ylabel('Magnitude (dB)', 'FontSize', 11);
title('Open-Loop Bode Plot -- Controller Comparison', 'FontSize', 13);
legend('Location','northeast','FontSize',9); grid on; hold off;
subplot(2,1,2);
yline(-180, 'k--', 'LineWidth', 1, 'Label', '-180 deg');
ylabel('Phase (deg)', 'FontSize', 11);
xlabel('Frequency (rad/s)', 'FontSize', 11);
legend('Location','southwest','FontSize',9); grid on; hold off;
saveas(gcf, 'figures/fig06_bode_controllers.png');
fprintf('Saved fig06_bode_controllers.png\n');

%% ── FIG 07: DISTURBANCE REJECTION ─────────────────────────────────────────
% Disturbance TF: from plant input to output = G/(1+C*G)
t_dist = 0:0.1:300;
d_signal = zeros(size(t_dist));
d_signal(t_dist >= 100) = 0.25;   % 25% step disturbance at t=100s

figure('Name','Disturbance Rejection','Position',[100 100 900 480]);
hold on;
for k = [2, 4, 5, 7]   % PI, PID, Lead, Lead-Lag
    L = controllers{k}*G;
    CL_r = feedback(L, 1);          % ref to output
    CL_d = feedback(G, controllers{k}); % disturbance to output
    [yr, ~] = step(CL_r, t_dist);
    [yd, ~] = lsim(CL_d, d_signal, t_dist);
    plot(t_dist, yr + yd, 'Color', colors{k}, 'LineStyle', styles{k}, ...
         'LineWidth', 1.8, 'DisplayName', ctrl_names{k});
end
xline(100, 'k:', 'LineWidth', 1.5, 'Label', 'Disturbance');
yline(1.0, 'k--', 'LineWidth', 1);
hold off;
xlabel('Time (s)', 'FontSize', 12); ylabel('Flow Rate (normalised)', 'FontSize', 12);
title('Disturbance Rejection -- 25% Step Disturbance at t=100 s', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 10);
ylim([0.5 1.6]); grid on; box on;
saveas(gcf, 'figures/fig07_disturbance.png');
fprintf('Saved fig07_disturbance.png\n');

%% ── FIG 08: PARAMETER SENSITIVITY (+/-20% T2) ───────────────────────────────
figure('Name','Parameter Sensitivity','Position',[100 100 1100 480]);
deltas = [0.8, 1.0, 1.2];
delta_labels = {'-20%', 'Nominal', '+20%'};
line_styles = {'--', '-', '-.'};

for ci = 1:2
    subplot(1,2,ci);
    hold on;
    ctrl_idx = [2, 4]; % PI, PID
    C_sel = controllers{ctrl_idx(ci)};
    for di = 1:3
        T2v = 7.05 * deltas(di);
        p2v = -1/T2v;
        num_v = [abs(p1)*abs(p2v)];
        den_v = [1, -(p1+p2v), p1*p2v];
        Gv = tf(num_v, den_v);
        CLv = feedback(C_sel*Gv, 1);
        [y,t] = step(CLv, t_sim);
        plot(t, y, 'LineStyle', line_styles{di}, 'LineWidth', 1.8, ...
             'DisplayName', ['T2 ' delta_labels{di}]);
    end
    yline(1.0, 'k:', 'LineWidth', 1);
    hold off;
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Flow Rate (normalised)', 'FontSize', 11);
    title([ctrl_names{ctrl_idx(ci)} ' Controller -- Parameter Sensitivity'], 'FontSize', 12);
    legend('FontSize', 9); grid on; xlim([0 200]);
end
saveas(gcf, 'figures/fig08_sensitivity.png');
fprintf('Saved fig08_sensitivity.png\n');

%% ── FIG 09: SENSOR NOISE ───────────────────────────────────────────────────
rng(42);
noise_std = 0.015;
noise = noise_std * randn(1, length(t_sim));

figure('Name','Sensor Noise','Position',[100 100 1100 480]);
for ci = 1:2
    subplot(1,2,ci);
    ctrl_idx = [2, 4];
    CL_sel = CL{ctrl_idx(ci)};
    [y_clean, ~] = lsim(CL_sel, ones(length(t_sim),1), t_sim);
    u_noisy = ones(length(t_sim),1) - noise';
    [y_noisy, ~] = lsim(CL_sel, u_noisy, t_sim);
    hold on;
    plot(t_sim, y_clean, 'b', 'LineWidth', 1.8, 'DisplayName', 'No noise');
    plot(t_sim, y_noisy, 'r', 'LineWidth', 0.9, 'DisplayName', ...
         sprintf('Noise \\sigma=%.1f%%', noise_std*100));
    hold off;
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Flow Rate (normalised)', 'FontSize', 11);
    title([ctrl_names{ctrl_idx(ci)} ' -- Sensor Noise Effect'], 'FontSize', 12);
    legend('FontSize', 9); grid on; xlim([0 200]);
end
saveas(gcf, 'figures/fig09_noise.png');
fprintf('Saved fig09_noise.png\n');

%% ── FIG 10: ACTUATOR SATURATION (Simulink-style using lsim approximation) ──
% Simulate saturation effect: when input is limited to u_max,
% the effective drive is reduced. We model this by scaling the
% closed-loop input signal and comparing responses.
figure('Name','Actuator Saturation','Position',[100 100 1100 480]);
u_max = 1.8;
for ci = 1:2
    subplot(1,2,ci);
    ctrl_idx = [2, 4];
    CL_sel = CL{ctrl_idx(ci)};
    % Normal response
    [y_lin, ~] = step(CL_sel, t_sim);
    % Approximation of saturated response: reduced reference input
    % (represents situation where actuator demand exceeds capacity)
    r_sat = min(ones(length(t_sim),1), u_max*0.6*ones(length(t_sim),1));
    [y_sat, ~] = lsim(CL_sel, r_sat, t_sim);
    hold on;
    plot(t_sim, y_lin, 'b', 'LineWidth', 2, 'DisplayName', 'Normal operation');
    plot(t_sim, y_sat, 'r--', 'LineWidth', 2, 'DisplayName', 'Saturation effect');
    yline(1.0, 'k:', 'LineWidth', 1);
    hold off;
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('Flow Rate (normalised)', 'FontSize', 11);
    title([ctrl_names{ctrl_idx(ci)} ' -- Actuator Saturation Effect'], 'FontSize', 12);
    legend('FontSize', 9); grid on; xlim([0 200]);
end
sgtitle('Actuator Saturation Effects (u_{max} = 1.8)', 'FontSize', 13);
saveas(gcf, 'figures/fig10_saturation.png');
fprintf('Saved fig10_saturation.png\n');

fprintf('\nScript 02 complete. All classical controller figures saved.\n');
