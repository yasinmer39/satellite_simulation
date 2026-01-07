% Attitude Estimation Test Script
% Faz 3: Yönelim Kestirimi - TRIAD, QUEST, MEKF test ve doğrulama
%
% Bu script, attitude estimation algoritmalarını test eder:
%   1. TRIAD algoritması testi
%   2. QUEST algoritması testi
%   3. MEKF simülasyonu (tam orbit)
%   4. Güneş konumu hesaplama testi
%   5. Kestirim vs Gerçek karşılaştırma (istatistiksel)
%
% Test Senaryosu:
%   - 60 kg LEO uydu, 1000 km yörünge
%   - IMU (STIM377H), Magnetometer (MAG-3), GNSS (OEM719)
%   - Sun sensor YOK - güneş konumu GNSS zamanından hesaplanır

clear; clc; close all;

fprintf('=== FAZ 3: YÖNELİM KESİTİRİMİ TEST ===\n\n');

%% ==================== SIMULATION PARAMETERS ====================
fprintf('1. Simülasyon Parametreleri...\n');

% Simülasyon süresi
t_sim = 600;         % 10 dakika (saniye)
dt = 0.1;            % 10 Hz (100 ms)
t = 0:dt:t_sim;
N = length(t);

% Başlangıç zamanı (UTC)
utc_start = datetime(2026, 6, 21, 12, 0, 0);  % Yaz gün dönümü

% Yörünge parametreleri (1000 km altitude, circular)
mu = 3.986004418e14;  % m^3/s^2
R_earth = 6378137;    % m
altitude = 1000e3;    % m
a = R_earth + altitude;
n = sqrt(mu / a^3);   % Mean motion (rad/s)
inc = 60 * pi/180;    % Inclination (rad)

%% ==================== TRUE ATTITUDE TRAJECTORY ====================
fprintf('2. Gerçek Yönelim Trajektorisi Oluşturuluyor...\n');

% Spacecraft angular velocity (yavaş rotasyon)
omega_body_true = [0.01; 0.005; 0.02] * pi/180;  % rad/s (çok yavaş spin)

% Initial attitude (Euler angles -> quaternion)
euler_init = [10; -5; 30] * pi/180;  % roll, pitch, yaw (rad)
q_true = zeros(4, N);
q_true(:,1) = euler_to_quat_local(euler_init);

% Propagate true attitude
for k = 2:N
    q_true(:,k) = propagate_quaternion(q_true(:,k-1), omega_body_true, dt);
end

fprintf('   - Başlangıç Euler açıları: [%.1f, %.1f, %.1f] deg\n', ...
        euler_init*180/pi);
fprintf('   - Angular velocity: [%.3f, %.3f, %.3f] deg/s\n', ...
        omega_body_true*180/pi);

%% ==================== ORBIT PROPAGATION ====================
fprintf('3. Yörünge Propagasyonu...\n');

r_eci = zeros(3, N);  % Satellite position (ECI)
v_eci = zeros(3, N);  % Satellite velocity (ECI)

% Initial position (circular orbit)
r_eci(:,1) = [a; 0; 0];
v_eci(:,1) = [0; sqrt(mu/a)*cos(inc); sqrt(mu/a)*sin(inc)];

% Simple Keplerian propagation
for k = 2:N
    mean_anomaly = n * t(k);
    
    % Position in orbital plane
    r_orbital = a * [cos(mean_anomaly); sin(mean_anomaly); 0];
    
    % Rotate to ECI (simplified: only inclination)
    R_inc = [1, 0, 0;
             0, cos(inc), -sin(inc);
             0, sin(inc), cos(inc)];
    r_eci(:,k) = R_inc * r_orbital;
    
    % Velocity
    v_orbital = a * n * [-sin(mean_anomaly); cos(mean_anomaly); 0];
    v_eci(:,k) = R_inc * v_orbital;
end

%% ==================== ENVIRONMENT MODELS ====================
fprintf('4. Çevre Modelleri (Manyetik Alan, Güneş)...\n');

% Reference vectors in ECI
B_eci = zeros(3, N);      % Magnetic field (ECI)
s_eci = zeros(3, N);      % Sun vector (ECI)
g_eci = zeros(3, N);      % Gravity (nadir) vector (ECI)

for k = 1:N
    % Current UTC time
    utc_current = utc_start + seconds(t(k));
    
    % Sun position (from almanac - deterministic)
    [s_eci(:,k), ~, ~, ~] = sun_position_eci_simple(utc_current);
    
    % Magnetic field (IGRF simplified - dipole model)
    B_eci(:,k) = igrf_dipole_simple(r_eci(:,k));
    
    % Gravity (nadir) direction
    g_eci(:,k) = -r_eci(:,k) / norm(r_eci(:,k));
end

%% ==================== SENSOR MEASUREMENTS (with noise) ====================
fprintf('5. Sensör Ölçümleri Simülasyonu...\n');

% Sensor noise parameters (from STIM377H, MAG-3)
gyro_noise_std = 0.15 * (pi/180) / sqrt(3600);  % ARW: 0.15 deg/sqrt(h)
gyro_bias_true = [0.5; -0.3; 0.4] * (pi/180) / 3600;  % 0.5 deg/h bias
accel_noise_std = 100e-6 * 9.81;  % 100 μg
mag_noise_std = 15e-9 / 50e-6;  % 15 nT / 50 μT = ~0.3% of full scale

% Measurement arrays
gyro_meas = zeros(3, N);
accel_meas = zeros(3, N);
mag_meas = zeros(3, N);

for k = 1:N
    % True attitude DCM (ECI -> Body)
    A_true = quat_to_dcm_local(q_true(:,k));
    
    % Gyro measurement (with bias and noise)
    gyro_meas(:,k) = omega_body_true + gyro_bias_true + ...
                     gyro_noise_std * sqrt(dt) * randn(3,1);
    
    % Accelerometer measurement (gravity in body frame)
    g_body_true = A_true * g_eci(:,k);
    accel_meas(:,k) = -g_body_true * 9.81 + accel_noise_std * randn(3,1);
    
    % Magnetometer measurement (B field in body frame)
    B_body_true = A_true * B_eci(:,k);
    mag_meas(:,k) = B_body_true + mag_noise_std * norm(B_body_true) * randn(3,1);
end

%% ==================== TEST 1: TRIAD ALGORITHM ====================
fprintf('\n6. Test 1: TRIAD Algoritması...\n');

q_triad = zeros(4, N);
triad_errors = zeros(1, N);

for k = 1:N
    % Normalize measurements
    accel_norm = -accel_meas(:,k) / norm(accel_meas(:,k));  % Points to nadir
    mag_norm = mag_meas(:,k) / norm(mag_meas(:,k));
    
    % Reference vectors
    g_ref = g_eci(:,k);
    B_ref = B_eci(:,k) / norm(B_eci(:,k));
    
    % TRIAD (gravity more accurate than mag)
    [q_triad(:,k), ~] = triad_algorithm(accel_norm, mag_norm, g_ref, B_ref);
    
    % Error calculation
    triad_errors(k) = quat_angle_error_local(q_triad(:,k), q_true(:,k));
end

fprintf('   - TRIAD Mean Error: %.3f deg\n', mean(triad_errors));
fprintf('   - TRIAD Std Error:  %.3f deg\n', std(triad_errors));
fprintf('   - TRIAD Max Error:  %.3f deg\n', max(triad_errors));

%% ==================== TEST 2: QUEST ALGORITHM ====================
fprintf('\n7. Test 2: QUEST Algoritması...\n');

q_quest = zeros(4, N);
quest_errors = zeros(1, N);

for k = 1:N
    % Measurements (body frame)
    accel_norm = -accel_meas(:,k) / norm(accel_meas(:,k));
    mag_norm = mag_meas(:,k) / norm(mag_meas(:,k));
    
    v_body = [accel_norm, mag_norm];
    
    % Reference vectors (ECI)
    g_ref = g_eci(:,k);
    B_ref = B_eci(:,k) / norm(B_eci(:,k));
    
    v_ref = [g_ref, B_ref];
    
    % Weights (gravity more accurate)
    weights = [0.7, 0.3];
    
    % QUEST
    [q_quest(:,k), ~, ~] = quest_algorithm(v_body, v_ref, weights);
    
    % Error calculation
    quest_errors(k) = quat_angle_error_local(q_quest(:,k), q_true(:,k));
end

fprintf('   - QUEST Mean Error: %.3f deg\n', mean(quest_errors));
fprintf('   - QUEST Std Error:  %.3f deg\n', std(quest_errors));
fprintf('   - QUEST Max Error:  %.3f deg\n', max(quest_errors));

%% ==================== TEST 3: MEKF ====================
fprintf('\n8. Test 3: MEKF (Multiplicative EKF)...\n');

q_mekf = zeros(4, N);
mekf_errors = zeros(1, N);
bias_est = zeros(3, N);

% Initialize MEKF
mekf_state = [];
P_mekf = [];

% MEKF parameters
mekf_params = struct();
mekf_params.gyro_noise = gyro_noise_std;
mekf_params.gyro_bias_noise = 0.5 * (pi/180) / 3600;
mekf_params.accel_noise = 0.03;
mekf_params.mag_noise = 0.05;
mekf_params.init_attitude_var = (10 * pi/180)^2;
mekf_params.init_bias_var = (1 * pi/180 / 3600)^2;

for k = 1:N
    % Measurements
    accel_norm = -accel_meas(:,k) / norm(accel_meas(:,k));
    mag_norm = mag_meas(:,k) / norm(mag_meas(:,k));
    
    % Reference vectors
    B_ref = B_eci(:,k) / norm(B_eci(:,k));
    g_ref = g_eci(:,k);
    
    % MEKF update
    [q_mekf(:,k), mekf_state, P_mekf] = mekf_attitude(...
        gyro_meas(:,k), accel_norm, mag_norm, ...
        B_ref, g_ref, dt, mekf_state, P_mekf, mekf_params);
    
    % Store bias estimate
    bias_est(:,k) = mekf_state.gyro_bias;
    
    % Error calculation
    mekf_errors(k) = quat_angle_error_local(q_mekf(:,k), q_true(:,k));
end

fprintf('   - MEKF Mean Error: %.3f deg\n', mean(mekf_errors));
fprintf('   - MEKF Std Error:  %.3f deg\n', std(mekf_errors));
fprintf('   - MEKF Max Error:  %.3f deg\n', max(mekf_errors));

% Bias estimation accuracy
bias_error = vecnorm(bias_est - gyro_bias_true) * 3600 * 180/pi;  % deg/h
fprintf('   - Bias Est Error (final): %.4f deg/h\n', bias_error(end));

%% ==================== TEST 4: SUN POSITION ====================
fprintf('\n9. Test 4: Güneş Konumu Hesaplama...\n');

% Calculate sun position at different times
test_times = [datetime(2026,1,1,12,0,0);   % Winter
              datetime(2026,3,21,12,0,0);  % Equinox
              datetime(2026,6,21,12,0,0);  % Summer
              datetime(2026,9,23,12,0,0)]; % Equinox

fprintf('   Güneş pozisyonları (ECI unit vector):\n');
for i = 1:length(test_times)
    [s, ~, ~, dist] = sun_position_eci_simple(test_times(i));
    fprintf('   %s: [%.4f, %.4f, %.4f], dist=%.4f AU\n', ...
            datestr(test_times(i), 'yyyy-mm-dd'), s(1), s(2), s(3), dist);
end

%% ==================== TEST 5: SUN VECTOR IN BODY FRAME ====================
fprintf('\n10. Test 5: Güneş Vektörü Body Frame (Kestirilen)...\n');

s_body_est = zeros(3, N);
s_body_true = zeros(3, N);
sun_angle_errors = zeros(1, N);

for k = 1:N
    % True sun vector in body frame
    A_true = quat_to_dcm_local(q_true(:,k));
    s_body_true(:,k) = A_true * s_eci(:,k);
    
    % Estimated sun vector in body frame (using MEKF quaternion)
    A_est = quat_to_dcm_local(q_mekf(:,k));
    s_body_est(:,k) = A_est * s_eci(:,k);
    
    % Angular error
    cos_angle = dot(s_body_true(:,k), s_body_est(:,k)) / ...
                (norm(s_body_true(:,k)) * norm(s_body_est(:,k)));
    sun_angle_errors(k) = acosd(min(1, max(-1, cos_angle)));
end

fprintf('   - Sun Vector Mean Error: %.3f deg\n', mean(sun_angle_errors));
fprintf('   - Sun Vector Std Error:  %.3f deg\n', std(sun_angle_errors));
fprintf('   - Sun Vector Max Error:  %.3f deg\n', max(sun_angle_errors));

%% ==================== STATISTICAL ANALYSIS ====================
fprintf('\n11. İstatistiksel Analiz (Gerçek vs Kestirim)...\n');

fprintf('\n   === KARŞILAŞTIRMA TABLOSU ===\n');
fprintf('   %-12s %10s %10s %10s\n', 'Algoritma', 'Mean(deg)', 'Std(deg)', 'Max(deg)');
fprintf('   %-12s %10.4f %10.4f %10.4f\n', 'TRIAD', mean(triad_errors), std(triad_errors), max(triad_errors));
fprintf('   %-12s %10.4f %10.4f %10.4f\n', 'QUEST', mean(quest_errors), std(quest_errors), max(quest_errors));
fprintf('   %-12s %10.4f %10.4f %10.4f\n', 'MEKF', mean(mekf_errors), std(mekf_errors), max(mekf_errors));

% 3-sigma bounds
fprintf('\n   3-Sigma Bounds:\n');
fprintf('   - TRIAD: %.3f deg\n', mean(triad_errors) + 3*std(triad_errors));
fprintf('   - QUEST: %.3f deg\n', mean(quest_errors) + 3*std(quest_errors));
fprintf('   - MEKF:  %.3f deg\n', mean(mekf_errors) + 3*std(mekf_errors));

%% ==================== PLOTS ====================
fprintf('\n12. Grafikler Oluşturuluyor...\n');

figure('Name', 'Attitude Estimation Results', 'Position', [100, 100, 1200, 800]);

% Plot 1: Attitude Errors
subplot(2,2,1);
plot(t, triad_errors, 'b-', 'LineWidth', 0.5); hold on;
plot(t, quest_errors, 'r-', 'LineWidth', 0.5);
plot(t, mekf_errors, 'g-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Attitude Error (deg)');
title('Attitude Estimation Error Comparison');
legend('TRIAD', 'QUEST', 'MEKF', 'Location', 'best');
grid on;

% Plot 2: MEKF Error Detail
subplot(2,2,2);
plot(t, mekf_errors, 'g-', 'LineWidth', 1); hold on;
plot(t, mean(mekf_errors)*ones(size(t)), 'k--', 'LineWidth', 1);
plot(t, (mean(mekf_errors)+3*std(mekf_errors))*ones(size(t)), 'r--', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('MEKF Error (deg)');
title('MEKF Attitude Error');
legend('Error', 'Mean', '3\sigma', 'Location', 'best');
grid on;

% Plot 3: Gyro Bias Estimation
subplot(2,2,3);
bias_true_dph = gyro_bias_true * 3600 * 180/pi;  % deg/h
plot(t, bias_est(1,:)*3600*180/pi, 'r-', 'LineWidth', 1); hold on;
plot(t, bias_est(2,:)*3600*180/pi, 'g-', 'LineWidth', 1);
plot(t, bias_est(3,:)*3600*180/pi, 'b-', 'LineWidth', 1);
plot([t(1), t(end)], [bias_true_dph(1), bias_true_dph(1)], 'r--');
plot([t(1), t(end)], [bias_true_dph(2), bias_true_dph(2)], 'g--');
plot([t(1), t(end)], [bias_true_dph(3), bias_true_dph(3)], 'b--');
xlabel('Time (s)');
ylabel('Gyro Bias (deg/h)');
title('Gyro Bias Estimation');
legend('b_x est', 'b_y est', 'b_z est', 'b_x true', 'b_y true', 'b_z true');
grid on;

% Plot 4: Sun Vector Error
subplot(2,2,4);
plot(t, sun_angle_errors, 'm-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Sun Vector Error (deg)');
title('Sun Direction Estimation Error (Body Frame)');
grid on;

saveas(gcf, 'attitude_estimation_results.png');
fprintf('   Grafik kaydedildi: attitude_estimation_results.png\n');

% Error histogram
figure('Name', 'Error Distribution', 'Position', [100, 100, 800, 400]);

subplot(1,2,1);
histogram(mekf_errors, 30, 'Normalization', 'probability');
xlabel('MEKF Error (deg)');
ylabel('Probability');
title('MEKF Error Distribution');
grid on;

subplot(1,2,2);
histogram(sun_angle_errors, 30, 'Normalization', 'probability');
xlabel('Sun Vector Error (deg)');
ylabel('Probability');
title('Sun Direction Error Distribution');
grid on;

saveas(gcf, 'error_distribution.png');
fprintf('   Grafik kaydedildi: error_distribution.png\n');

fprintf('\n=== TEST TAMAMLANDI ===\n');

%% ==================== LOCAL HELPER FUNCTIONS ====================

function q = euler_to_quat_local(euler)
    roll = euler(1); pitch = euler(2); yaw = euler(3);
    cr = cos(roll/2); sr = sin(roll/2);
    cp = cos(pitch/2); sp = sin(pitch/2);
    cy = cos(yaw/2); sy = sin(yaw/2);
    q = [sr*cp*cy - cr*sp*sy;
         cr*sp*cy + sr*cp*sy;
         cr*cp*sy - sr*sp*cy;
         cr*cp*cy + sr*sp*sy];
    q = q / norm(q);
end

function q_new = propagate_quaternion(q, omega, dt)
    wx = omega(1); wy = omega(2); wz = omega(3);
    Omega = [0, wz, -wy, wx;
            -wz, 0, wx, wy;
            wy, -wx, 0, wz;
            -wx, -wy, -wz, 0];
    q_dot = 0.5 * Omega * q;
    q_new = q + q_dot * dt;
    q_new = q_new / norm(q_new);
end

function A = quat_to_dcm_local(q)
    q1=q(1); q2=q(2); q3=q(3); q4=q(4);
    A = [q1^2-q2^2-q3^2+q4^2, 2*(q1*q2+q3*q4), 2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4), -q1^2+q2^2-q3^2+q4^2, 2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4), 2*(q2*q3-q1*q4), -q1^2-q2^2+q3^2+q4^2];
end

function angle = quat_angle_error_local(q_est, q_true)
    q_err = quat_mult_local(q_true, quat_inv_local(q_est));
    if q_err(4) < 0, q_err = -q_err; end
    angle = 2 * acos(min(1, abs(q_err(4)))) * 180/pi;
end

function q_inv = quat_inv_local(q)
    q_inv = [-q(1); -q(2); -q(3); q(4)];
end

function r = quat_mult_local(p, q)
    p1=p(1); p2=p(2); p3=p(3); p4=p(4);
    q1=q(1); q2=q(2); q3=q(3); q4=q(4);
    r = [p4*q1+p1*q4+p2*q3-p3*q2;
         p4*q2-p1*q3+p2*q4+p3*q1;
         p4*q3+p1*q2-p2*q1+p3*q4;
         p4*q4-p1*q1-p2*q2-p3*q3];
end

function [s_eci, s_km, eclipse, dist] = sun_position_eci_simple(utc_time)
    if isa(utc_time, 'datetime')
        year = utc_time.Year; month = utc_time.Month; day = utc_time.Day;
        hour = utc_time.Hour; minute = utc_time.Minute; second = utc_time.Second;
    else
        year = utc_time(1); month = utc_time(2); day = utc_time(3);
        hour = utc_time(4); minute = utc_time(5); second = utc_time(6);
    end
    day_frac = day + (hour + minute/60 + second/3600) / 24;
    JD = 367*year - floor(7*(year + floor((month+9)/12))/4) + floor(275*month/9) + day_frac + 1721013.5;
    T = (JD - 2451545.0) / 36525;
    
    lambda_M = mod(280.4606184 + 36000.77005361*T, 360);
    M = mod(357.5277233 + 35999.05034*T, 360);
    M_rad = M * pi/180;
    lambda = mod(lambda_M + 1.914666471*sin(M_rad) + 0.019994643*sin(2*M_rad), 360);
    lambda_rad = lambda * pi/180;
    dist = 1.000140612 - 0.016708617*cos(M_rad) - 0.000139589*cos(2*M_rad);
    epsilon = (23.439291 - 0.0130042*T) * pi/180;
    
    s_eci = [cos(lambda_rad); cos(epsilon)*sin(lambda_rad); sin(epsilon)*sin(lambda_rad)];
    s_eci = s_eci / norm(s_eci);
    s_km = s_eci * dist * 149597870.7;
    eclipse = 0;
end

function B_eci = igrf_dipole_simple(r_eci)
    % Simple dipole model for testing
    B0 = 3.12e-5;  % Tesla at equator surface
    R_earth = 6378137;  % m
    r = norm(r_eci);
    r_hat = r_eci / r;
    
    % Dipole moment (aligned with z-axis for simplicity)
    m_hat = [0; 0; 1];
    
    % Dipole field
    B_eci = B0 * (R_earth/r)^3 * (3 * dot(m_hat, r_hat) * r_hat - m_hat);
end
