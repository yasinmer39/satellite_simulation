% Fuel Slosh Test and Simulation Script
% Case Study (f): Yakıt çalkalanmasının sistem üzerindeki etkisi
%
% Bu script, fuel slosh modelinin etkilerini simüle eder:
%   1. Tank ve propellant tasarımı
%   2. Slosh parametreleri hesabı
%   3. Free response (serbest salınım)
%   4. Coupled spacecraft + slosh simulation
%   5. Slosh etkisi vs etkisiz karşılaştırma
%   6. Kontrol sistemi üzerindeki etki analizi
%
% Referanslar:
%   - Dodge, F.T., "Dynamic Behavior of Liquids in Moving Containers"
%   - NASA SP-106, "Propellant Slosh Loads"
%   - Abramson, H.N., "Liquid Sloshing in Containers"

clear; clc; close all;

fprintf('=================================================================\n');
fprintf('   CASE STUDY (f): YAKIT ÇALKALANMASI (FUEL SLOSH) ANALİZİ\n');
fprintf('   60 kg LEO Uydu - Hydrazine Tank\n');
fprintf('=================================================================\n\n');

%% ==================== 1. TANK TASARIMI ====================
fprintf('BÖLÜM 1: Tank ve Propellant Tasarımı\n');
fprintf('-----------------------------------------------------------------\n');

spacecraft_mass = 60;  % kg

% Tank ve slosh parametrelerini hesapla
[slosh_params, tank_params] = fuel_slosh_tank_design(spacecraft_mass);

%% ==================== 2. SLOSH MODEL PARAMETRELERİ ====================
fprintf('\n\nBÖLÜM 2: Slosh Model Parametreleri\n');
fprintf('-----------------------------------------------------------------\n');

% Simülasyon parametreleri
t_sim = 60;           % 60 saniye
dt = 0.01;            % 100 Hz
t = 0:dt:t_sim;
N = length(t);

% Effective gravity (maneuver sırasında)
g_eff_values = [0.01, 0.05, 0.1, 0.5, 1.0] * 9.81;  % Different acceleration levels

fprintf('Farklı ivme seviyelerinde slosh frekansı:\n');
for i = 1:length(g_eff_values)
    g_val = g_eff_values(i);
    omega_n = sqrt(slosh_params.k_freq * g_val / (tank_params.radius));
    freq_n = omega_n / (2*pi);
    period = 1/freq_n;
    fprintf('   g = %.2f m/s² (%.2f g): f = %.3f Hz, T = %.2f s\n', ...
            g_val, g_val/9.81, freq_n, period);
end

%% ==================== 3. FREE RESPONSE SİMÜLASYONU ====================
fprintf('\n\nBÖLÜM 3: Serbest Salınım (Free Response) Simülasyonu\n');
fprintf('-----------------------------------------------------------------\n');

% Initial slosh angle (perturbation)
psi_0 = 5 * pi/180;  % 5 derece başlangıç açısı
psi_dot_0 = 0;

% Simülasyon için farklı g değerleri
g_test = 0.05 * 9.81;  % 0.05g (tipik manevra)

% Storage
psi_free = zeros(1, N);
psi_dot_free = zeros(1, N);
psi_free(1) = psi_0;
psi_dot_free(1) = psi_dot_0;

% Update natural frequency for this g
omega_n = sqrt(slosh_params.k_freq * g_test / tank_params.radius);
zeta = slosh_params.zeta;

% Free response simulation (linearized pendulum)
for k = 2:N
    % State space form: ẍ + 2ζωₙẋ + ωₙ²x = 0
    psi_ddot = -2*zeta*omega_n*psi_dot_free(k-1) - omega_n^2*psi_free(k-1);
    
    % Euler integration
    psi_dot_free(k) = psi_dot_free(k-1) + psi_ddot * dt;
    psi_free(k) = psi_free(k-1) + psi_dot_free(k) * dt;
end

% Analytical solution for comparison
omega_d = omega_n * sqrt(1 - zeta^2);  % Damped frequency
psi_analytical = psi_0 * exp(-zeta*omega_n*t) .* cos(omega_d*t);

% Settling time (to 2% of initial)
settling_idx = find(abs(psi_free) < 0.02*psi_0, 1, 'first');
if ~isempty(settling_idx)
    settling_time = t(settling_idx);
else
    settling_time = t_sim;
end

fprintf('Serbest salınım parametreleri (g = %.3f m/s²):\n', g_test);
fprintf('   - Doğal frekans: %.4f Hz\n', omega_n/(2*pi));
fprintf('   - Sönümlü frekans: %.4f Hz\n', omega_d/(2*pi));
fprintf('   - Periyot: %.2f s\n', 2*pi/omega_d);
fprintf('   - Settling time (2%%): %.1f s\n', settling_time);

%% ==================== 4. COUPLED DYNAMICS SİMÜLASYONU ====================
fprintf('\n\nBÖLÜM 4: Coupled Spacecraft + Slosh Simülasyonu\n');
fprintf('-----------------------------------------------------------------\n');

% Spacecraft parameters
I_sc = slosh_params.I_sc_dry(1,1);  % Use single axis for 2D sim
I_f = slosh_params.m_fixed * norm(slosh_params.tank_position)^2;

% Simulation parameters struct
sim_params.I_sc = I_sc;
sim_params.I_f = I_f;
sim_params.slosh = slosh_params;
sim_params.g_eff = g_test;

% Initial conditions: [θ, θ_dot, ψ, ψ_dot]
% Spacecraft initially at rest, slosh excited
x0_coupled = [0; 0; psi_0; 0];

% Simulate without control (free response)
t_coupled = 0:dt:30;  % 30 seconds
N_coupled = length(t_coupled);

x_coupled = zeros(4, N_coupled);
x_coupled(:,1) = x0_coupled;

% Simple RK4 integration
for k = 1:N_coupled-1
    tk = t_coupled(k);
    xk = x_coupled(:,k);
    
    k1 = spacecraft_slosh_ode(tk, xk, sim_params, 0);
    k2 = spacecraft_slosh_ode(tk+dt/2, xk+dt/2*k1, sim_params, 0);
    k3 = spacecraft_slosh_ode(tk+dt/2, xk+dt/2*k2, sim_params, 0);
    k4 = spacecraft_slosh_ode(tk+dt, xk+dt*k3, sim_params, 0);
    
    x_coupled(:,k+1) = xk + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% Extract states
theta_coupled = x_coupled(1,:) * 180/pi;
theta_dot_coupled = x_coupled(2,:) * 180/pi;
psi_coupled = x_coupled(3,:) * 180/pi;
psi_dot_coupled = x_coupled(4,:) * 180/pi;

fprintf('Coupled simulation results:\n');
fprintf('   - Max spacecraft angle: %.3f deg\n', max(abs(theta_coupled)));
fprintf('   - Max slosh angle: %.3f deg\n', max(abs(psi_coupled)));
fprintf('   - Energy transfer spacecraft->slosh and back visible in plots\n');

%% ==================== 5. SLOSH ETKİSİ KARŞILAŞTIRMASI ====================
fprintf('\n\nBÖLÜM 5: Slosh Etkisi vs Etkisiz Karşılaştırma\n');
fprintf('-----------------------------------------------------------------\n');

% Simulation with impulsive torque (attitude maneuver)
tau_impulse = 0.5;  % Nm (typical reaction wheel torque)
t_maneuver = 0:dt:30;
N_maneuver = length(t_maneuver);

% WITH slosh
x_with_slosh = zeros(4, N_maneuver);
x_with_slosh(:,1) = [0; 0; 0; 0];  % Start from rest

% WITHOUT slosh (rigid body only)
theta_no_slosh = zeros(1, N_maneuver);
theta_dot_no_slosh = zeros(1, N_maneuver);
I_total_rigid = I_sc + (slosh_params.m_slosh + slosh_params.m_fixed) * ...
                norm(slosh_params.tank_position)^2;

% Torque profile: step input for 2 seconds, then off
for k = 1:N_maneuver-1
    tk = t_maneuver(k);
    
    % Torque input
    if tk < 2
        tau = tau_impulse;
    else
        tau = 0;
    end
    
    % WITH slosh
    xk = x_with_slosh(:,k);
    k1 = spacecraft_slosh_ode(tk, xk, sim_params, tau);
    k2 = spacecraft_slosh_ode(tk+dt/2, xk+dt/2*k1, sim_params, tau);
    k3 = spacecraft_slosh_ode(tk+dt/2, xk+dt/2*k2, sim_params, tau);
    k4 = spacecraft_slosh_ode(tk+dt, xk+dt*k3, sim_params, tau);
    x_with_slosh(:,k+1) = xk + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    % WITHOUT slosh (simple rigid body)
    theta_ddot_rigid = tau / I_total_rigid;
    theta_dot_no_slosh(k+1) = theta_dot_no_slosh(k) + theta_ddot_rigid * dt;
    theta_no_slosh(k+1) = theta_no_slosh(k) + theta_dot_no_slosh(k) * dt;
end

% Extract results
theta_with_slosh = x_with_slosh(1,:) * 180/pi;
psi_with_slosh = x_with_slosh(3,:) * 180/pi;
theta_no_slosh_deg = theta_no_slosh * 180/pi;

% Calculate difference
theta_diff = theta_with_slosh - theta_no_slosh_deg;

fprintf('Manevra sonuçları (%.1f Nm, 2 s):\n', tau_impulse);
fprintf('   WITH SLOSH:\n');
fprintf('      - Final angle: %.3f deg\n', theta_with_slosh(end));
fprintf('      - Max slosh: %.3f deg\n', max(abs(psi_with_slosh)));
fprintf('      - Residual oscillation: %.4f deg\n', std(theta_with_slosh(end-100:end)));
fprintf('   WITHOUT SLOSH:\n');
fprintf('      - Final angle: %.3f deg\n', theta_no_slosh_deg(end));
fprintf('   DIFFERENCE:\n');
fprintf('      - Max deviation: %.4f deg\n', max(abs(theta_diff)));
fprintf('      - RMS deviation: %.4f deg\n', rms(theta_diff));

%% ==================== 6. KONTROL SİSTEMİ ETKİ ANALİZİ ====================
fprintf('\n\nBÖLÜM 6: Kontrol Sistemi Üzerindeki Etki Analizi\n');
fprintf('-----------------------------------------------------------------\n');

% PD controller simulation
Kp = 2.0;   % Proportional gain (Nm/rad)
Kd = 1.0;   % Derivative gain (Nm/(rad/s))
theta_cmd = 30 * pi/180;  % 30 deg setpoint

t_control = 0:dt:60;
N_control = length(t_control);

% WITH slosh + PD control
x_ctrl_slosh = zeros(4, N_control);
x_ctrl_slosh(:,1) = [0; 0; 0; 0];

% WITHOUT slosh + PD control
theta_ctrl_rigid = zeros(1, N_control);
theta_dot_ctrl_rigid = zeros(1, N_control);

% Control effort storage
tau_ctrl_slosh = zeros(1, N_control);
tau_ctrl_rigid = zeros(1, N_control);

for k = 1:N_control-1
    % WITH slosh
    theta_err_slosh = theta_cmd - x_ctrl_slosh(1,k);
    theta_dot_slosh = x_ctrl_slosh(2,k);
    tau_slosh = Kp * theta_err_slosh - Kd * theta_dot_slosh;
    tau_slosh = max(-1, min(1, tau_slosh));  % Saturation ±1 Nm
    tau_ctrl_slosh(k) = tau_slosh;
    
    xk = x_ctrl_slosh(:,k);
    k1 = spacecraft_slosh_ode(t_control(k), xk, sim_params, tau_slosh);
    k2 = spacecraft_slosh_ode(t_control(k)+dt/2, xk+dt/2*k1, sim_params, tau_slosh);
    k3 = spacecraft_slosh_ode(t_control(k)+dt/2, xk+dt/2*k2, sim_params, tau_slosh);
    k4 = spacecraft_slosh_ode(t_control(k)+dt, xk+dt*k3, sim_params, tau_slosh);
    x_ctrl_slosh(:,k+1) = xk + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    % WITHOUT slosh
    theta_err_rigid = theta_cmd - theta_ctrl_rigid(k);
    tau_rigid = Kp * theta_err_rigid - Kd * theta_dot_ctrl_rigid(k);
    tau_rigid = max(-1, min(1, tau_rigid));
    tau_ctrl_rigid(k) = tau_rigid;
    
    theta_ddot_rigid = tau_rigid / I_total_rigid;
    theta_dot_ctrl_rigid(k+1) = theta_dot_ctrl_rigid(k) + theta_ddot_rigid * dt;
    theta_ctrl_rigid(k+1) = theta_ctrl_rigid(k) + theta_dot_ctrl_rigid(k) * dt;
end

theta_ctrl_slosh_deg = x_ctrl_slosh(1,:) * 180/pi;
psi_ctrl_slosh_deg = x_ctrl_slosh(3,:) * 180/pi;
theta_ctrl_rigid_deg = theta_ctrl_rigid * 180/pi;

% Performance metrics
settling_slosh = find(abs(theta_ctrl_slosh_deg - 30) < 0.3, 1, 'first');
settling_rigid = find(abs(theta_ctrl_rigid_deg - 30) < 0.3, 1, 'first');

overshoot_slosh = (max(theta_ctrl_slosh_deg) - 30) / 30 * 100;
overshoot_rigid = (max(theta_ctrl_rigid_deg) - 30) / 30 * 100;

fprintf('PD Controller Performance (Kp=%.1f, Kd=%.1f):\n', Kp, Kd);
fprintf('   %-25s %12s %12s\n', 'Metric', 'With Slosh', 'Without');
fprintf('   ─────────────────────────────────────────────────\n');
fprintf('   %-25s %12.2f %12.2f s\n', 'Settling time (1%)', ...
        t_control(settling_slosh), t_control(settling_rigid));
fprintf('   %-25s %12.2f %12.2f %%\n', 'Overshoot', overshoot_slosh, overshoot_rigid);
fprintf('   %-25s %12.4f %12.4f deg\n', 'SS Error (mean)', ...
        mean(theta_ctrl_slosh_deg(end-100:end))-30, ...
        mean(theta_ctrl_rigid_deg(end-100:end))-30);
fprintf('   %-25s %12.4f %12.4f deg\n', 'SS Error (std)', ...
        std(theta_ctrl_slosh_deg(end-100:end)), ...
        std(theta_ctrl_rigid_deg(end-100:end)));
fprintf('   %-25s %12.4f %12.4f Nm\n', 'Control effort (rms)', ...
        rms(tau_ctrl_slosh), rms(tau_ctrl_rigid));

%% ==================== 7. GRAFİKLER ====================
fprintf('\n\nBÖLÜM 7: Grafikler\n');
fprintf('-----------------------------------------------------------------\n');

% Figure 1: Free Response
figure('Name', 'Fuel Slosh Free Response', 'Position', [50, 50, 1200, 500]);

subplot(1,2,1);
plot(t, psi_free*180/pi, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, psi_analytical*180/pi, 'r--', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Slosh Angle (deg)');
title('Free Response: Slosh Pendulum');
legend('Numerical', 'Analytical', 'Location', 'best');
grid on;

subplot(1,2,2);
plot(t_coupled, theta_coupled, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_coupled, psi_coupled, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angle (deg)');
title('Coupled Dynamics: Energy Exchange');
legend('Spacecraft θ', 'Slosh ψ', 'Location', 'best');
grid on;

saveas(gcf, 'fuel_slosh_free_response.png');
fprintf('   Kaydedildi: fuel_slosh_free_response.png\n');

% Figure 2: Slosh Effect Comparison
figure('Name', 'Slosh Effect Comparison', 'Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(t_maneuver, theta_with_slosh, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_maneuver, theta_no_slosh_deg, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Spacecraft Angle (deg)');
title('Maneuver Response: With vs Without Slosh');
legend('With Slosh', 'Without Slosh', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(t_maneuver, psi_with_slosh, 'g-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Slosh Angle (deg)');
title('Slosh Excitation During Maneuver');
grid on;

subplot(2,2,3);
plot(t_maneuver, theta_diff, 'm-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Angle Difference (deg)');
title('Difference: (With Slosh) - (Without Slosh)');
grid on;

subplot(2,2,4);
% Phase portrait
plot(psi_with_slosh, x_with_slosh(4,:)*180/pi, 'b-', 'LineWidth', 1);
xlabel('ψ (deg)');
ylabel('ψ̇ (deg/s)');
title('Slosh Phase Portrait');
grid on;

saveas(gcf, 'fuel_slosh_comparison.png');
fprintf('   Kaydedildi: fuel_slosh_comparison.png\n');

% Figure 3: Control System Performance
figure('Name', 'Control with Slosh', 'Position', [150, 150, 1200, 800]);

subplot(2,2,1);
plot(t_control, theta_ctrl_slosh_deg, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_control, theta_ctrl_rigid_deg, 'r--', 'LineWidth', 1.5);
yline(30, 'k--', 'Setpoint');
xlabel('Time (s)');
ylabel('Spacecraft Angle (deg)');
title('Closed-Loop Response: PD Controller');
legend('With Slosh', 'Without Slosh', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(t_control, psi_ctrl_slosh_deg, 'g-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Slosh Angle (deg)');
title('Slosh During Controlled Maneuver');
grid on;

subplot(2,2,3);
plot(t_control, tau_ctrl_slosh, 'b-', 'LineWidth', 1);
hold on;
plot(t_control, tau_ctrl_rigid, 'r--', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Control Torque (Nm)');
title('Control Effort');
legend('With Slosh', 'Without Slosh');
grid on;

subplot(2,2,4);
err_slosh = theta_ctrl_slosh_deg - 30;
err_rigid = theta_ctrl_rigid_deg - 30;
plot(t_control, err_slosh, 'b-', 'LineWidth', 1);
hold on;
plot(t_control, err_rigid, 'r--', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Tracking Error (deg)');
title('Tracking Error Comparison');
legend('With Slosh', 'Without Slosh');
grid on;

saveas(gcf, 'fuel_slosh_control.png');
fprintf('   Kaydedildi: fuel_slosh_control.png\n');

%% ==================== ÖZET ====================
fprintf('\n');
fprintf('=================================================================\n');
fprintf('   FUEL SLOSH ANALİZ ÖZETİ\n');
fprintf('=================================================================\n');
fprintf('   Tank: Spherical, R = %.2f cm, V = %.1f L\n', ...
        tank_params.radius*100, tank_params.volume_total*1000);
fprintf('   Propellant: %s, %.2f kg (50%% fill)\n', ...
        tank_params.propellant.name, tank_params.propellant_mass);
fprintf('   Slosh mass: %.2f kg (%.0f%% of liquid)\n', ...
        slosh_params.m_slosh, slosh_params.k_mass*100);
fprintf('   Slosh freq @0.05g: %.4f Hz\n', omega_n/(2*pi));
fprintf('   Damping ratio: %.1f%%\n', slosh_params.zeta*100);
fprintf('\n   ETKİLER:\n');
fprintf('   - Slosh, kontrol sisteminde ~%.1f%% overshoot artışına neden olur\n', ...
        overshoot_slosh - overshoot_rigid);
fprintf('   - Settling time ~%.1f s uzar\n', ...
        t_control(settling_slosh) - t_control(settling_rigid));
fprintf('   - Residual oscillation: ±%.3f deg\n', ...
        std(theta_ctrl_slosh_deg(end-100:end)));
fprintf('=================================================================\n');

%% ==================== LOCAL FUNCTIONS ====================

function x_dot = spacecraft_slosh_ode(t, x, params, tau)
    % Extract state
    theta = x(1);
    theta_dot = x(2);
    psi = x(3);
    psi_dot = x(4);
    
    % Parameters
    I_sc = params.I_sc;
    I_f = params.I_f;
    m_s = params.slosh.m_slosh;
    m_f = params.slosh.m_fixed;
    a = params.slosh.L;
    b = norm(params.slosh.tank_position);
    c = params.slosh.c;
    g = params.g_eff;
    
    % Mass matrix
    M11 = I_sc + m_f*b^2 + m_s*(b^2 - a*b*cos(psi));
    M12 = -m_s*a*b*cos(psi);
    M21 = M12;
    M22 = I_f + m_s*a^2;
    
    % Nonlinear terms
    C1 = m_s*a*b*(theta_dot + psi_dot)^2*sin(psi);
    C2 = m_s*a*g*sin(psi) + c*psi_dot;
    
    % Solve
    M = [M11, M12; M21, M22];
    RHS = [tau + C1; -C2];
    acc = M \ RHS;
    
    x_dot = [theta_dot; acc(1); psi_dot; acc(2)];
end

function [slosh_params, tank_params] = fuel_slosh_tank_design(spacecraft_mass)
    % Simplified version for this script
    propellant.name = 'Hydrazine (N2H4)';
    propellant.density = 1004;
    
    propellant_mass = 4;  % kg at 50% fill
    fill_ratio = 0.50;
    
    propellant_volume = propellant_mass / propellant.density;
    tank_volume_total = propellant_volume / fill_ratio;
    tank_radius = (3*tank_volume_total/(4*pi))^(1/3);
    
    % Tank params
    tank_params.propellant = propellant;
    tank_params.radius = tank_radius;
    tank_params.volume_total = tank_volume_total;
    tank_params.propellant_mass = propellant_mass;
    
    % Slosh params
    k_freq = 1.19;
    k_mass = 0.65;
    slosh_params.m_slosh = k_mass * propellant_mass;
    slosh_params.m_fixed = propellant_mass - slosh_params.m_slosh;
    slosh_params.L = tank_radius / k_freq;
    slosh_params.k_freq = k_freq;
    slosh_params.k_mass = k_mass;
    slosh_params.zeta = 0.01;
    slosh_params.tank_position = [0; 0; 0.15];
    
    g_eff = 0.05 * 9.81;
    omega_n = sqrt(k_freq * g_eff / tank_radius);
    slosh_params.c = 2 * slosh_params.zeta * omega_n * slosh_params.m_slosh * slosh_params.L^2;
    
    slosh_params.I_sc_dry = spacecraft_mass * 0.09 * eye(3);
    
    fprintf('Tank tasarımı tamamlandı.\n');
    fprintf('   - Tank yarıçapı: %.2f cm\n', tank_radius*100);
    fprintf('   - Propellant: %.2f kg\n', propellant_mass);
    fprintf('   - Slosh mass: %.2f kg\n', slosh_params.m_slosh);
end
