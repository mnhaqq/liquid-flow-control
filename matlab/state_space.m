%% CPEN 401 -- Group 16: Liquid Flow Control in a Pipeline Pumping Station
%% Script 03: State-Space Design -- Pole Placement, Observer, LQR
%  Generates:
%    fig11_state_feedback.png
%    fig12_observer.png

clear; clc; close all;

%% ── PLANT ──────────────────────────────────────────────────────────────────
T1 = 2.0; T2 = 7.05;
p1 = -1/T1; p2 = -1/T2;
num_G = [abs(p1)*abs(p2)];
den_G = [1, -(p1+p2), p1*p2];
G = tf(num_G, den_G);
t_sim = 0:0.05:200;

% Reference PID for comparison
C_PID = tf([3.0, 8.0, 0.4], [1, 0]);
CL_PID = feedback(C_PID*G, 1);
[y_PID, ~] = step(CL_PID, t_sim);

%% ── STATE-SPACE MODEL ──────────────────────────────────────────────────────
% States: x1 = Pp (pump pressure, normalised), x2 = Q (flow)
A = [p1,   0;
     1.0,  p2];
B = [abs(p1); 0];
C_out = [0, 1];
D_out = 0;

sys_ss = ss(A, B, C_out, D_out);
fprintf('=== STATE-SPACE MODEL ===\n');
fprintf('A = \n'); disp(A);
fprintf('B = \n'); disp(B);
fprintf('C = [%g %g]\n\n', C_out);

%% ── CONTROLLABILITY & OBSERVABILITY ───────────────────────────────────────
Mc = ctrb(A, B);
Mo = obsv(A, C_out);
rank_c = rank(Mc);
rank_o = rank(Mo);
fprintf('Controllability matrix rank = %d/2 -> %s\n', rank_c, ...
        ternary(rank_c==2,'FULLY CONTROLLABLE','NOT controllable'));
fprintf('Observability matrix rank   = %d/2 -> %s\n', rank_o, ...
        ternary(rank_o==2,'FULLY OBSERVABLE','NOT observable'));

%% ── POLE PLACEMENT ─────────────────────────────────────────────────────────
des_poles_1 = [-0.6, -1.0];   % moderate
des_poles_2 = [-1.2, -2.0];   % aggressive

K_pp1 = place(A, B, des_poles_1);
K_pp2 = place(A, B, des_poles_2);
fprintf('\n=== POLE PLACEMENT ===\n');
fprintf('Desired poles set 1: [%.1f, %.1f]\n', des_poles_1);
fprintf('Feedback gain K1 = [%.5f, %.5f]\n', K_pp1);
fprintf('Desired poles set 2: [%.1f, %.1f]\n', des_poles_2);
fprintf('Feedback gain K2 = [%.5f, %.5f]\n', K_pp2);

% Pre-compensator: ensures unity DC gain
A_cl1 = A - B*K_pp1;  N1 = -1 / (C_out*(A_cl1\B));
A_cl2 = A - B*K_pp2;  N2 = -1 / (C_out*(A_cl2\B));
fprintf('Precompensator N1 = %.4f,  N2 = %.4f\n', N1, N2);

% Simulate state feedback using ss with augmented input
sys_pp1 = ss(A_cl1, B*N1, C_out, 0);
sys_pp2 = ss(A_cl2, B*N2, C_out, 0);
[y_pp1, ~] = step(sys_pp1, t_sim);
[y_pp2, ~] = step(sys_pp2, t_sim);

% Print metrics
for k = 1:2
    if k==1, y=y_pp1; label='PP Set 1'; else, y=y_pp2; label='PP Set 2'; end
    yss = y(end);
    idx10 = find(y>=0.10*yss,1); idx90 = find(y>=0.90*yss,1);
    if ~isempty(idx10)&&~isempty(idx90), tr=t_sim(idx90)-t_sim(idx10); else, tr=NaN; end
    idx_out = find(abs(y-yss)>0.02*yss); ts = t_sim(idx_out(end));
    fprintf('%s: tr=%.2fs  ts=%.2fs  OS=%.1f%%\n', label, tr, ts, ...
            max(0,(max(y)-yss)/yss*100));
end

%% ── LQR ────────────────────────────────────────────────────────────────────
Q_lqr = diag([1.0, 50.0]);
R_lqr = 0.05;
[K_lqr, ~, eig_lqr] = lqr(A, B, Q_lqr, R_lqr);
fprintf('\n=== LQR ===\n');
fprintf('K_lqr = [%.5f, %.5f]\n', K_lqr);
fprintf('Closed-loop eigenvalues: %.4f, %.4f\n', eig_lqr(1), eig_lqr(2));

A_lqr = A - B*K_lqr;
N_lqr = -1 / (C_out*(A_lqr\B));
sys_lqr = ss(A_lqr, B*N_lqr, C_out, 0);
[y_lqr, ~] = step(sys_lqr, t_sim);

%% ── FIG 11: STATE FEEDBACK COMPARISON ─────────────────────────────────────
figure('Name','State Feedback','Position',[100 100 900 520]);
hold on;
plot(t_sim, y_pp1, 'b',  'LineWidth', 2.0, ...
     'DisplayName', sprintf('Pole Placement [%.1f, %.1f]', des_poles_1));
plot(t_sim, y_pp2, 'r',  'LineWidth', 2.0, ...
     'DisplayName', sprintf('Pole Placement [%.1f, %.1f]', des_poles_2));
plot(t_sim, y_lqr, 'g--','LineWidth', 2.0, 'DisplayName', 'LQR');
plot(t_sim, y_PID, 'm:', 'LineWidth', 1.8, 'DisplayName', 'PID (classical ref.)');
yline(1.0, 'k:', 'LineWidth', 1.2);
hold off;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Flow Rate (normalised)', 'FontSize', 12);
title('State Feedback -- Pole Placement & LQR vs PID', 'FontSize', 13);
legend('Location','southeast','FontSize',10);
xlim([0 200]); grid on; box on;
saveas(gcf, 'figures/fig11_state_feedback.png');
fprintf('\nSaved fig11_state_feedback.png\n');

%% ── LUENBERGER OBSERVER ─────────────────────────────────────────────────────
obs_poles = [-3.0, -4.5];
L_obs = place(A', C_out', obs_poles)';
fprintf('\n=== OBSERVER ===\n');
fprintf('Observer poles: [%.1f, %.1f]\n', obs_poles);
fprintf('Observer gain L = [%.5f; %.5f]\n', L_obs);

% Simulate observer-based feedback
% Combined system: [x_true; x_hat] with u = -K*x_hat + N*r
A_obs_true = A - B*K_pp1;          % closed with true state (ideal)
A_obs_hat  = A - L_obs*C_out;      % observer dynamics

% Simulate using ode45 for the coupled system
r_in = 1.0;
x0 = zeros(4,1);  % [x_true; x_hat]

[t_obs, X] = ode45(@(t,x) obs_dynamics(t, x, A, B, C_out, K_pp1, N1, L_obs, r_in), ...
                   [0, 200], x0);

y_true_obs = X(:,2);           % x2 = Q (true flow)
obs_error  = vecnorm(X(:,1:2) - X(:,3:4), 2, 2);  % ||x - x_hat||

%% ── FIG 12: OBSERVER DESIGN ────────────────────────────────────────────────
figure('Name','Observer','Position',[100 100 900 680]);
subplot(2,1,1);
hold on;
plot(t_sim, y_pp1, 'b--', 'LineWidth', 1.6, 'DisplayName', 'Full state feedback');
plot(t_obs, y_true_obs, 'r', 'LineWidth', 1.8, 'DisplayName', 'Observer-based feedback');
yline(1.0, 'k:', 'LineWidth', 1);
hold off;
xlabel('Time (s)', 'FontSize', 11);
ylabel('Flow Rate (normalised)', 'FontSize', 11);
title('Observer-Based State Feedback', 'FontSize', 12);
legend('FontSize', 10); grid on; xlim([0 200]);

subplot(2,1,2);
semilogy(t_obs, obs_error + 1e-12, 'g', 'LineWidth', 1.8);
xlabel('Time (s)', 'FontSize', 11);
ylabel('Observer Error |x - xhat| (log)', 'FontSize', 11);
title('State Estimation Error Convergence', 'FontSize', 12);
grid on; xlim([0 200]);
saveas(gcf, 'figures/fig12_observer.png');
fprintf('Saved fig12_observer.png\n');

fprintf('\nScript 03 complete.\n');

%% ── HELPER FUNCTIONS ───────────────────────────────────────────────────────
function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function dxdt = obs_dynamics(~, x, A, B, C, K, N, L, r)
    x_true = x(1:2);
    x_hat  = x(3:4);
    u = -K*x_hat + N*r;
    y_meas = C*x_true;
    y_hat  = C*x_hat;
    dx_true = A*x_true + B*u;
    dx_hat  = A*x_hat  + B*u + L*(y_meas - y_hat);
    dxdt = [dx_true; dx_hat];
end
