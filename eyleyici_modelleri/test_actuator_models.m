% =========================================================================
% Eyleyici Modelleri Test Scripti
% 60 kg LEO Uydu
% =========================================================================
% Bu script, tüm eyleyici modellerini test eder:
%   1. Reaction Wheel (tek tekerlek)
%   2. Reaction Wheel Assembly (3 eksen)
%   3. Magnetorquer (tek eksen)
%   4. Magnetorquer Assembly (3 eksen + B-dot)
%   5. Gas Jet (thruster)
% =========================================================================

clear; clc; close all;

fprintf('========================================\n');
fprintf('   EYLEYİCİ MODELLERİ TESTİ\n');
fprintf('   60 kg LEO Uydu\n');
fprintf('========================================\n\n');

%% Simülasyon Parametreleri
dt = 0.01;          % Zaman adımı (s) - 100 Hz
t_end = 60;         % Simülasyon süresi (s)
t = 0:dt:t_end;
N = length(t);

fprintf('Simülasyon parametreleri:\n');
fprintf('  dt = %.3f s, T = %.1f s, N = %d samples\n\n', dt, t_end, N);

%% TEST 1: Tek Reaction Wheel
fprintf('TEST 1: TEK REACTION WHEEL\n');
fprintf('---------------------------\n');

wheel_state = [];
speed_history = zeros(1, N);
torque_history = zeros(1, N);
h_history = zeros(1, N);

for i = 1:N
    % Sinüsoidal torque komutu
    torque_cmd = 0.01 * sin(2*pi*0.05*t(i));  % 0.01 N·m, 0.05 Hz
    
    [wheel_state, torque_actual, ~] = reaction_wheel_model(wheel_state, torque_cmd, dt, []);
    
    speed_history(i) = wheel_state.speed;
    torque_history(i) = torque_actual;
    h_history(i) = wheel_state.h;
end

fprintf('Reaction Wheel sonuçları:\n');
fprintf('  Max hız: %.1f rpm\n', max(abs(speed_history)));
fprintf('  Max momentum: %.4f N·m·s\n', max(abs(h_history)));
fprintf('  Max tork: %.4f N·m\n\n', max(abs(torque_history)));

%% TEST 2: 3-Axis Reaction Wheel Assembly
fprintf('TEST 2: 3-AXIS REACTION WHEEL ASSEMBLY\n');
fprintf('---------------------------------------\n');

rwa_state = [];
omega_sc = [0.01; 0.01; 0.05];  % rad/s başlangıç açısal hızı
I_sc = diag([5, 4, 3]);        % kg·m² spacecraft inertia

omega_history = zeros(3, N);
h_rwa_history = zeros(3, N);
torque_rwa_history = zeros(3, N);

for i = 1:N
    % Basit PD kontrol
    omega_desired = [0; 0; 0];
    Kp = 0.1;
    Kd = 0.5;
    
    torque_cmd = -Kp * omega_sc - Kd * omega_sc;  % Basitleştirilmiş
    
    [rwa_state, torque_actual, h_wheels, ~] = reaction_wheel_assembly(rwa_state, torque_cmd, dt, []);
    
    % Spacecraft dynamics (simplified)
    omega_dot = I_sc \ (torque_actual - cross(omega_sc, I_sc * omega_sc));
    omega_sc = omega_sc + omega_dot * dt;
    
    omega_history(:,i) = omega_sc;
    h_rwa_history(:,i) = h_wheels;
    torque_rwa_history(:,i) = torque_actual;
end

fprintf('RWA sonuçları:\n');
fprintf('  Final ω: [%.4f, %.4f, %.4f] rad/s\n', omega_sc(1), omega_sc(2), omega_sc(3));
fprintf('  Final h_rwa: [%.4f, %.4f, %.4f] N·m·s\n', h_rwa_history(1,end), h_rwa_history(2,end), h_rwa_history(3,end));
fprintf('  Damping achieved: %.1f%%\n\n', (1 - norm(omega_sc)/norm([0.01;0.01;0.05]))*100);

%% TEST 3: Tek Magnetorquer
fprintf('TEST 3: TEK MAGNETORQUER\n');
fprintf('-------------------------\n');

% Sabit manyetik alan (body frame)
B_body = [20000; 10000; 40000] * 1e-9;  % T (20, 10, 40 µT)

dipole_cmd = [3; 0; 0];  % A·m² (sadece X ekseni)

[dipole_actual, torque_mtq, power_mtq] = magnetorquer_model(dipole_cmd, B_body, []);

fprintf('Magnetorquer sonuçları:\n');
fprintf('  Komut dipol: [%.2f, %.2f, %.2f] A·m²\n', dipole_cmd(1), dipole_cmd(2), dipole_cmd(3));
fprintf('  Gerçek dipol: [%.3f, %.3f, %.3f] A·m²\n', dipole_actual(1), dipole_actual(2), dipole_actual(3));
fprintf('  Tork: [%.2e, %.2e, %.2e] N·m\n', torque_mtq(1), torque_mtq(2), torque_mtq(3));
fprintf('  Güç: %.2f W\n\n', power_mtq);

%% TEST 4: B-dot Detumbling with Magnetorquer Assembly
fprintf('TEST 4: B-DOT DETUMBLING\n');
fprintf('-------------------------\n');

omega_sc_bdot = [0.1; 0.1; 0.1];  % rad/s (yüksek başlangıç hızı)
I_sc_bdot = diag([5, 4, 3]);

omega_bdot_history = zeros(3, N);
B_body_prev = [20000; 10000; 40000] * 1e-9;

for i = 1:N
    % Yörünge boyunca değişen manyetik alan simülasyonu
    orbit_phase = 2*pi * t(i) / 5400;  % ~90 dakika yörünge
    B_inertial = [30000 * cos(orbit_phase); 
                  20000 * sin(orbit_phase); 
                  40000] * 1e-9;  % T
    
    % Basit DCM (sadece z ekseni etrafında dönüş)
    theta_z = omega_sc_bdot(3) * t(i);
    Rz = [cos(theta_z), sin(theta_z), 0;
         -sin(theta_z), cos(theta_z), 0;
          0, 0, 1];
    B_body_current = Rz * B_inertial;
    
    % B-dot control
    [dipole_actual, torque_mtq, ~] = magnetorquer_assembly([], B_body_current, [], 'bdot', B_body_prev, dt);
    
    % Spacecraft dynamics
    omega_dot = I_sc_bdot \ (torque_mtq - cross(omega_sc_bdot, I_sc_bdot * omega_sc_bdot));
    omega_sc_bdot = omega_sc_bdot + omega_dot * dt;
    
    omega_bdot_history(:,i) = omega_sc_bdot;
    B_body_prev = B_body_current;
end

fprintf('B-dot detumbling sonuçları:\n');
fprintf('  Başlangıç |ω|: %.4f rad/s\n', norm([0.1; 0.1; 0.1]));
fprintf('  Final |ω|: %.4f rad/s\n', norm(omega_sc_bdot));
fprintf('  Azalma oranı: %.1f%%\n\n', (1 - norm(omega_sc_bdot)/norm([0.1;0.1;0.1]))*100);

%% TEST 5: Gas Jet Thruster
fprintf('TEST 5: GAS JET THRUSTER\n');
fprintf('-------------------------\n');

jet_state = [];
thrust_history = zeros(1, N);
firing_history = false(1, N);

% Pulse modulated firing pattern
pulse_period = 2.0;  % s
pulse_duty = 0.3;    % 30% duty cycle

for i = 1:N
    % PWM pattern
    phase = mod(t(i), pulse_period) / pulse_period;
    cmd_on = phase < pulse_duty;
    
    [jet_state, thrust, ~] = gas_jet_model(cmd_on, t(i), jet_state, []);
    
    thrust_history(i) = thrust;
    firing_history(i) = jet_state.firing;
end

fprintf('Gas Jet sonuçları:\n');
fprintf('  Nominal thrust: %.2f N\n', 1.0);
fprintf('  Total impulse: %.3f N·s\n', sum(thrust_history) * dt);
fprintf('  Total firing time: %.2f s\n', sum(firing_history) * dt);
fprintf('  Effective duty cycle: %.1f%%\n\n', sum(firing_history)/N*100);

%% Grafikler
fprintf('GRAFİKLER OLUŞTURULUYOR...\n');

figure('Name', 'Eyleyici Modelleri Test Sonuçları', 'Position', [50 50 1400 900]);

% 1. Reaction Wheel Speed
subplot(3,4,1);
plot(t, speed_history, 'b-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('Hız (rpm)');
title('Reaction Wheel Hızı');
grid on;

% 2. Reaction Wheel Torque
subplot(3,4,2);
plot(t, torque_history * 1000, 'r-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('Tork (mN·m)');
title('Reaction Wheel Torku');
grid on;

% 3. Reaction Wheel Momentum
subplot(3,4,3);
plot(t, h_history * 1000, 'g-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('Momentum (mN·m·s)');
title('Reaction Wheel Momentumu');
grid on;

% 4. RWA - Spacecraft Angular Velocity
subplot(3,4,4);
plot(t, omega_history(1,:)*180/pi, 'r-', 'LineWidth', 1);
hold on;
plot(t, omega_history(2,:)*180/pi, 'g-', 'LineWidth', 1);
plot(t, omega_history(3,:)*180/pi, 'b-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('ω (°/s)');
title('RWA ile Hız Sönümleme');
legend('ω_x', 'ω_y', 'ω_z', 'Location', 'best');
grid on;

% 5. RWA - Wheel Momentum
subplot(3,4,5);
plot(t, h_rwa_history(1,:)*1000, 'r-', 'LineWidth', 1);
hold on;
plot(t, h_rwa_history(2,:)*1000, 'g-', 'LineWidth', 1);
plot(t, h_rwa_history(3,:)*1000, 'b-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('h (mN·m·s)');
title('RWA Tekerlek Momentumları');
legend('h_x', 'h_y', 'h_z', 'Location', 'best');
grid on;

% 6. RWA - Applied Torque
subplot(3,4,6);
plot(t, torque_rwa_history(1,:)*1000, 'r-', 'LineWidth', 1);
hold on;
plot(t, torque_rwa_history(2,:)*1000, 'g-', 'LineWidth', 1);
plot(t, torque_rwa_history(3,:)*1000, 'b-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('Tork (mN·m)');
title('RWA Uygulanan Tork');
legend('τ_x', 'τ_y', 'τ_z', 'Location', 'best');
grid on;

% 7. B-dot Detumbling - Angular Velocity
subplot(3,4,7);
plot(t, omega_bdot_history(1,:)*180/pi, 'r-', 'LineWidth', 1);
hold on;
plot(t, omega_bdot_history(2,:)*180/pi, 'g-', 'LineWidth', 1);
plot(t, omega_bdot_history(3,:)*180/pi, 'b-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('ω (°/s)');
title('B-dot Detumbling');
legend('ω_x', 'ω_y', 'ω_z', 'Location', 'best');
grid on;

% 8. B-dot - Angular Velocity Magnitude
subplot(3,4,8);
omega_mag = vecnorm(omega_bdot_history);
semilogy(t, omega_mag*180/pi, 'b-', 'LineWidth', 1.5);
xlabel('Zaman (s)');
ylabel('|ω| (°/s)');
title('B-dot |ω| Azalması');
grid on;

% 9. Gas Jet - Thrust Profile
subplot(3,4,9);
plot(t, thrust_history, 'b-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('Thrust (N)');
title('Gas Jet Thrust Profili');
ylim([-0.1, 1.5]);
grid on;

% 10. Gas Jet - Zoom on transient
subplot(3,4,10);
idx = t < 5;  % İlk 5 saniye
plot(t(idx), thrust_history(idx), 'b-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('Thrust (N)');
title('Gas Jet Rise/Fall (Zoom)');
grid on;

% 11. Gas Jet - Cumulative Impulse
subplot(3,4,11);
cumulative_impulse = cumsum(thrust_history) * dt;
plot(t, cumulative_impulse, 'b-', 'LineWidth', 1.5);
xlabel('Zaman (s)');
ylabel('Kümülatif İmpuls (N·s)');
title('Gas Jet Toplam İmpuls');
grid on;

% 12. Özet Tablo
subplot(3,4,12);
axis off;
text(0.1, 0.95, 'EYLEYİCİ ÖZETİ', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.80, 'Reaction Wheel:', 'FontSize', 10, 'FontWeight', 'bold');
text(0.1, 0.70, sprintf('  Max tork: 20 mN·m'), 'FontSize', 9);
text(0.1, 0.60, sprintf('  Max hız: 6000 rpm'), 'FontSize', 9);
text(0.1, 0.50, 'Magnetorquer:', 'FontSize', 10, 'FontWeight', 'bold');
text(0.1, 0.40, sprintf('  Max dipol: 5 A·m²'), 'FontSize', 9);
text(0.1, 0.30, 'Gas Jet:', 'FontSize', 10, 'FontWeight', 'bold');
text(0.1, 0.20, sprintf('  Nominal: 1.0 N'), 'FontSize', 9);
text(0.1, 0.10, sprintf('  Isp: 70 s'), 'FontSize', 9);

sgtitle('Eyleyici Modelleri Test Sonuçları', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nTüm testler tamamlandı.\n');
fprintf('========================================\n');
