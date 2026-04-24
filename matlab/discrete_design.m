%% CPEN 401 -- Group 16: Liquid Flow Control in a Pipeline Pumping Station
%% Script 04: Discrete-Time Controller Design
%  Generates:
%    fig13_discrete.png
%    fig14_nyquist.png
%    fig15_states.png

clear; clc; close all;

%% ── PLANT ──────────────────────────────────────────────────────────────────
T1 = 2.0; T2 = 7.05;
p1 = -1/T1; p2 = -1/T2;
num_G = [abs(p1)*abs(p2)];
den_G = [1, -(p1+p2), p1*p2];
G = tf(num_G, den_G);

A = [p1, 0; 1.0, p2];
B = [abs(p1); 0];
C_out = [0, 1]; D_out = 0;
sys_ss = ss(A, B, C_out, D_out);

% Sampling period
Ts = 1.0;   % seconds

fprintf('=== DISCRETE-TIME DESIGN ===\n');
fprintf('Sampling period Ts = %.1f s\n', Ts);
fprintf('Nyquist frequency  = %.4f rad/s\n', pi/Ts);

%% ── DISCRETISE PLANT (ZOH) ─────────────────────────────────────────────────
sys_dt = c2d(sys_ss, Ts, 'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sys_dt);
eig_d = eig(Ad);

fprintf('\nDiscrete A matrix:\n'); disp(Ad);
fprintf('Discrete B matrix:\n'); disp(Bd);
fprintf('Discrete eigenvalues: %.6f, %.6f\n', eig_d(1), eig_d(2));
fprintf('All |lambda| < 1: %d (stable)\n', all(abs(eig_d) < 1));

%% ── DISCRETE POLE PLACEMENT ─────────────────────────────────────────────────
% Map continuous desired poles to z-domain: z = exp(lambda*Ts)
des_c1 = [-0.6, -1.0]; des_c2 = [-1.2, -2.0];
z_des1 = exp(des_c1 * Ts);
z_des2 = exp(des_c2 * Ts);
fprintf('\nDiscrete desired poles set 1: %.4f, %.4f\n', z_des1);
fprintf('Discrete desired poles set 2: %.4f, %.4f\n', z_des2);

K_d1 = place(Ad, Bd, z_des1);
K_d2 = place(Ad, Bd, z_des2);
fprintf('Gain K_d1 = [%.5f, %.5f]\n', K_d1);
fprintf('Gain K_d2 = [%.5f, %.5f]\n', K_d2);

% Pre-compensator for discrete system
Nd1 = 1 / (Cd * ((eye(2) - Ad + Bd*K_d1) \ Bd));
Nd2 = 1 / (Cd * ((eye(2) - Ad + Bd*K_d2) \ Bd));

% Simulate discrete pole placement
n_steps = 201;
t_d = (0:n_steps-1) * Ts;



y_d1 = sim_discrete_sf(Ad, Bd, Cd, K_d1, Nd1, n_steps);
y_d2 = sim_discrete_sf(Ad, Bd, Cd, K_d2, Nd2, n_steps);

%% ── DISCRETE PID (TUSTIN) ───────────────────────────────────────────────────
Kp=8.0; Ki=0.4; Kd_=7.0;
b0 = Kp + Ki*Ts/2 + 2*Kd_/Ts;
b1 = -Kp + Ki*Ts/2 - 4*Kd_/Ts;
b2 = 2*Kd_/Ts;
fprintf('\nDiscrete PID Tustin coefficients:\n');
fprintf('  b0 = %.4f,  b1 = %.4f,  b2 = %.4f\n', b0, b1, b2);

% Simulate discrete PID on continuous plant (using ode45 between steps)
y_dpid = zeros(n_steps,1);
x_plant = zeros(2,1);
e_prev1 = 0; e_prev2 = 0; u_prev = 0;
for k = 1:n_steps
    y_k = C_out * x_plant;
    y_dpid(k) = y_k;
    e_k = 1.0 - y_k;
    u_k = u_prev + b0*e_k + b1*e_prev1 + b2*e_prev2;
    u_k = max(-5, min(5, u_k));  % saturation
    % Propagate plant for Ts seconds
    [~, X_ode] = ode45(@(t,x) A*x + B*u_k, [0 Ts], x_plant);
    x_plant = X_ode(end,:)';
    e_prev2 = e_prev1; e_prev1 = e_k; u_prev = u_k;
end

%% ── CONTINUOUS REFERENCES ───────────────────────────────────────────────────
des_p1 = [-0.6, -1.0];
K_pp1 = place(A, B, des_p1);
A_cl1 = A - B*K_pp1;
N1 = -1/(C_out*(A_cl1\B));
sys_pp1 = ss(A_cl1, B*N1, C_out, 0);
t_cont = 0:0.1:200;
[y_pp1_cont, ~] = step(sys_pp1, t_cont);

C_PID = tf([3.0, 8.0, 0.4], [1, 0]);
CL_PID = feedback(C_PID*G, 1);
[y_PID_cont, ~] = step(CL_PID, t_cont);

%% ── FIG 13: CONTINUOUS VS DISCRETE ─────────────────────────────────────────
figure('Name','Discrete vs Continuous','Position',[100 100 1100 520]);
subplot(1,2,1);
hold on;
plot(t_cont, y_pp1_cont, 'b', 'LineWidth', 2, 'DisplayName', 'Continuous PP');
stairs(t_d, y_d1, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Discrete PP (T_s=%gs)', Ts));
yline(1.0, 'k:', 'LineWidth', 1);
hold off;
xlabel('Time (s)','FontSize',11); ylabel('Flow Rate (normalised)','FontSize',11);
title('Pole Placement -- Continuous vs Discrete','FontSize',12);
legend('FontSize',9); grid on; xlim([0 200]);

subplot(1,2,2);
hold on;
plot(t_cont, y_PID_cont, 'b', 'LineWidth', 2, 'DisplayName', 'Continuous PID');
stairs(t_d, y_dpid, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Discrete PID Tustin (T_s=%gs)',Ts));
yline(1.0,'k:','LineWidth',1);
hold off;
xlabel('Time (s)','FontSize',11); ylabel('Flow Rate (normalised)','FontSize',11);
title('PID -- Continuous vs Discrete','FontSize',12);
legend('FontSize',9); grid on; xlim([0 200]);
saveas(gcf, 'figures/fig13_discrete.png');
fprintf('Saved fig13_discrete.png\n');

%% ── FIG 14: NYQUIST PLOTS ───────────────────────────────────────────────────
C_PI  = tf(6.0*[20, 1], [20, 0]);

figure('Name','Nyquist','Position',[100 100 1100 520]);
subplot(1,2,1);
nyquist(C_PI*G);
title('Nyquist Plot -- PI Controller','FontSize',12);
grid on;
xlim([-5 5]); ylim([-5 5]);

subplot(1,2,2);
nyquist(C_PID*G);
title('Nyquist Plot -- PID Controller','FontSize',12);
grid on;
xlim([-5 5]); ylim([-5 5]);
saveas(gcf, 'figures/fig14_nyquist.png');
fprintf('Saved fig14_nyquist.png\n');

%% ── FIG 15: STATE TRAJECTORIES ─────────────────────────────────────────────
des_p2 = [-1.2, -2.0];
K_pp2 = place(A, B, des_p2);
A_cl2 = A - B*K_pp2;
N2 = -1/(C_out*(A_cl2\B));


t_st = 0:0.1:200;
X1 = sim_states(A, B, K_pp1, N1, t_st);
X2 = sim_states(A, B, K_pp2, N2, t_st);

figure('Name','State Trajectories','Position',[100 100 1100 500]);
subplot(1,2,1);
hold on;
plot(t_st, X1(:,1), 'b', 'LineWidth',1.8, 'DisplayName', 'x_1 -- Pump pressure');
plot(t_st, X1(:,2), 'r', 'LineWidth',1.8, 'DisplayName', 'x_2 -- Flow rate');
hold off;
xlabel('Time (s)','FontSize',11); ylabel('State value','FontSize',11);
title(sprintf('State Trajectories -- Poles [%.1f, %.1f]',des_c1),'FontSize',12);
legend('FontSize',9); grid on; xlim([0 200]);

subplot(1,2,2);
hold on;
plot(t_st, X2(:,1), 'b', 'LineWidth',1.8, 'DisplayName', 'x_1 -- Pump pressure');
plot(t_st, X2(:,2), 'r', 'LineWidth',1.8, 'DisplayName', 'x_2 -- Flow rate');
hold off;
xlabel('Time (s)','FontSize',11); ylabel('State value','FontSize',11);
title(sprintf('State Trajectories -- Poles [%.1f, %.1f]',des_c2),'FontSize',12);
legend('FontSize',9); grid on; xlim([0 200]);
saveas(gcf, 'figures/fig15_states.png');
fprintf('Saved fig15_states.png\n');
fprintf('\nScript 04 complete. All discrete and final figures saved.\n');

%% ── HELPER FUNCTION ─────────────────────────────────────────────────────────
function y_hist = sim_discrete_sf(Ad, Bd, Cd, K, N, n_steps)
    x = zeros(size(Ad,1),1);
    y_hist = zeros(n_steps,1);
    for k = 1:n_steps
        y_hist(k) = Cd * x;
        u = -K*x + N;
        x = Ad*x + Bd*u;
    end
end

function X = sim_states(A, B, K, N, t_vec)
    dt = t_vec(2)-t_vec(1);
    x = zeros(2,1);
    X = zeros(length(t_vec), 2);
    for k = 1:length(t_vec)
        X(k,:) = x';
        u = -K*x + N;
        x = x + dt*(A*x + B*u);
    end
end
