% =========================================================================
% Sensör Modelleri Test Scripti
% 60 kg LEO Uydu
% =========================================================================
% Bu script, tüm sensör modellerini test eder:
%   1. STIM377H Gyroscope
%   2. STIM377H Accelerometer
%   3. MAG-3 Magnetometer
%   4. GNSS-701 GPS Receiver
% =========================================================================

clear; clc; close all;

fprintf('========================================\n');
fprintf('   SENSÖR MODELLERİ TESTİ\n');
fprintf('   60 kg LEO Uydu\n');
fprintf('========================================\n\n');

%% Simülasyon Parametreleri
dt = 0.01;          % Zaman adımı (s) - 100 Hz
t_end = 60;         % Simülasyon süresi (s)
t = 0:dt:t_end;
N = length(t);

fprintf('Simülasyon parametreleri:\n');
fprintf('  dt = %.3f s, T = %.1f s, N = %d samples\n\n', dt, t_end, N);

%% TEST 1: STIM377H Gyroscope
fprintf('TEST 1: STIM377H GYROSCOPE\n');
fprintf('----------------------------\n');

% Gerçek açısal hız (sabit + sinüsoidal pertürbasyon)
omega_true_base = [0.01; 0.005; 0.1];  % rad/s (yaklaşık 6°/s yaw)
omega_true_all = zeros(3, N);
omega_meas_all = zeros(3, N);

gyro_state = [];
params_gyro = [];  % varsayılan parametreler

for i = 1:N
    % Küçük sinüsoidal pertürbasyon ekle
    omega_true = omega_true_base + 0.001 * sin(2*pi*0.1*t(i)) * [1; 1; 0];
    omega_true_all(:,i) = omega_true;
    
    [omega_meas, gyro_state] = gyro_model_stim377h(omega_true, dt, params_gyro, gyro_state);
    omega_meas_all(:,i) = omega_meas;
end

% Hata hesapla
gyro_error = omega_meas_all - omega_true_all;
gyro_error_deg_h = gyro_error * 180/pi * 3600;  % °/h

fprintf('Gyro ölçüm hataları (°/h):\n');
fprintf('  X-axis: mean = %.2f, std = %.2f\n', mean(gyro_error_deg_h(1,:)), std(gyro_error_deg_h(1,:)));
fprintf('  Y-axis: mean = %.2f, std = %.2f\n', mean(gyro_error_deg_h(2,:)), std(gyro_error_deg_h(2,:)));
fprintf('  Z-axis: mean = %.2f, std = %.2f\n\n', mean(gyro_error_deg_h(3,:)), std(gyro_error_deg_h(3,:)));

%% TEST 2: STIM377H Accelerometer
fprintf('TEST 2: STIM377H ACCELEROMETER\n');
fprintf('-------------------------------\n');

g0 = 9.80665;  % m/s²

% Gerçek ivme (mikro-yerçekimi + pertürbasyonlar)
accel_true_base = [0; 0; -1e-6] * g0;  % Çok küçük yerçekimi etkisi
accel_true_all = zeros(3, N);
accel_meas_all = zeros(3, N);

accel_state = [];
params_accel = [];

for i = 1:N
    % Küçük pertürbasyon ekle (drag, SRP vb. simülasyonu)
    accel_true = accel_true_base + 1e-7 * g0 * sin(2*pi*0.05*t(i)) * [1; 0; 0];
    accel_true_all(:,i) = accel_true;
    
    [accel_meas, accel_state] = accel_model_stim377h(accel_true, dt, params_accel, accel_state);
    accel_meas_all(:,i) = accel_meas;
end

% Hata hesapla
accel_error = accel_meas_all - accel_true_all;
accel_error_mg = accel_error / g0 * 1000;  % milli-g

fprintf('Accelerometer ölçüm hataları (mg):\n');
fprintf('  X-axis: mean = %.4f, std = %.4f\n', mean(accel_error_mg(1,:)), std(accel_error_mg(1,:)));
fprintf('  Y-axis: mean = %.4f, std = %.4f\n', mean(accel_error_mg(2,:)), std(accel_error_mg(2,:)));
fprintf('  Z-axis: mean = %.4f, std = %.4f\n\n', mean(accel_error_mg(3,:)), std(accel_error_mg(3,:)));

%% TEST 3: MAG-3 Magnetometer
fprintf('TEST 3: MAG-3 MAGNETOMETER\n');
fprintf('---------------------------\n');

% Gerçek manyetik alan (LEO'da tipik değerler)
% ~30,000-50,000 nT @ 1000 km altitude
B_true_base = [20000; 10000; 40000];  % nT
B_true_all = zeros(3, N);
B_meas_all = zeros(3, N);

mag_state = [];
params_mag = [];

for i = 1:N
    % Yörünge boyunca değişen manyetik alan simülasyonu
    orbit_phase = 2*pi * t(i) / 5400;  % ~90 dakika yörünge periyodu
    B_true = B_true_base + 5000 * [cos(orbit_phase); sin(orbit_phase); 0.5*cos(2*orbit_phase)];
    B_true_all(:,i) = B_true;
    
    [B_meas, mag_state] = magnetometer_model_mag3(B_true, params_mag, mag_state);
    B_meas_all(:,i) = B_meas;
end

% Hata hesapla
mag_error = B_meas_all - B_true_all;

fprintf('Magnetometer ölçüm hataları (nT):\n');
fprintf('  X-axis: mean = %.1f, std = %.1f\n', mean(mag_error(1,:)), std(mag_error(1,:)));
fprintf('  Y-axis: mean = %.1f, std = %.1f\n', mean(mag_error(2,:)), std(mag_error(2,:)));
fprintf('  Z-axis: mean = %.1f, std = %.1f\n', mean(mag_error(3,:)), std(mag_error(3,:)));
fprintf('  |B| relative error: %.3f %%\n\n', ...
    mean(abs(vecnorm(B_meas_all) - vecnorm(B_true_all)) ./ vecnorm(B_true_all)) * 100);

%% TEST 4: GNSS-701 GPS Receiver
fprintf('TEST 4: GNSS-701 GPS RECEIVER\n');
fprintf('------------------------------\n');

% Gerçek yörünge (basit dairesel LEO)
R_earth = 6378.137;  % km
altitude = 1000;     % km
r_orbit = R_earth + altitude;
omega_orbit = sqrt(398600 / r_orbit^3);  % rad/s

% 1 Hz örnekleme (GPS update rate)
dt_gps = 1.0;
t_gps = 0:dt_gps:t_end;
N_gps = length(t_gps);

r_true_all = zeros(3, N_gps);
v_true_all = zeros(3, N_gps);
r_meas_all = zeros(3, N_gps);
v_meas_all = zeros(3, N_gps);
valid_all = false(1, N_gps);

gnss_state = [];
params_gnss = [];

for i = 1:N_gps
    % Dairesel yörünge
    theta = omega_orbit * t_gps(i);
    r_true = r_orbit * [cos(theta); sin(theta); 0];
    v_true = r_orbit * omega_orbit * [-sin(theta); cos(theta); 0];
    
    r_true_all(:,i) = r_true;
    v_true_all(:,i) = v_true;
    
    [r_meas, v_meas, ~, valid, gnss_state] = gnss_model_701(r_true, v_true, t_gps(i), params_gnss, gnss_state);
    
    if valid
        r_meas_all(:,i) = r_meas;
        v_meas_all(:,i) = v_meas;
    else
        r_meas_all(:,i) = NaN(3,1);
        v_meas_all(:,i) = NaN(3,1);
    end
    valid_all(i) = valid;
end

% Geçerli ölçümler için hata hesapla
valid_idx = find(valid_all);
if ~isempty(valid_idx)
    pos_error = r_meas_all(:,valid_idx) - r_true_all(:,valid_idx);
    pos_error_m = pos_error * 1000;  % km -> m
    
    vel_error = v_meas_all(:,valid_idx) - v_true_all(:,valid_idx);
    vel_error_ms = vel_error * 1000;  % km/s -> m/s
    
    fprintf('GPS ölçüm hataları:\n');
    fprintf('  Position 3D RMS: %.2f m\n', rms(vecnorm(pos_error_m)));
    fprintf('  Velocity 3D RMS: %.4f m/s\n', rms(vecnorm(vel_error_ms)));
    fprintf('  Fix availability: %.1f %%\n', sum(valid_all)/N_gps*100);
    fprintf('  First fix at: %.1f s\n\n', t_gps(valid_idx(1)));
else
    fprintf('GPS fix alınamadı!\n\n');
end

%% Grafikler
fprintf('GRAFİKLER OLUŞTURULUYOR...\n');

figure('Name', 'Sensör Modelleri Test Sonuçları', 'Position', [50 50 1400 900]);

% 1. Gyro Ölçümleri
subplot(3,4,1);
plot(t, omega_true_all(3,:)*180/pi, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, omega_meas_all(3,:)*180/pi, 'r-', 'LineWidth', 0.5);
xlabel('Zaman (s)');
ylabel('ω_z (°/s)');
title('STIM377H Gyro - Z Ekseni');
legend('Gerçek', 'Ölçüm', 'Location', 'best');
grid on;

% 2. Gyro Hataları
subplot(3,4,2);
plot(t, gyro_error_deg_h(1,:), 'r-', 'LineWidth', 0.5);
hold on;
plot(t, gyro_error_deg_h(2,:), 'g-', 'LineWidth', 0.5);
plot(t, gyro_error_deg_h(3,:), 'b-', 'LineWidth', 0.5);
xlabel('Zaman (s)');
ylabel('Hata (°/h)');
title('Gyro Ölçüm Hataları');
legend('X', 'Y', 'Z', 'Location', 'best');
grid on;

% 3. Gyro Allan Variance (basit tahmin)
subplot(3,4,3);
% Kümülatif açı hatası
angle_error = cumsum(gyro_error(3,:)) * dt * 180/pi;  % derece
plot(t, angle_error, 'b-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('Kümülatif Açı Hatası (°)');
title('Gyro Drift Etkisi');
grid on;

% 4. Accelerometer Ölçümleri
subplot(3,4,4);
plot(t, accel_true_all(1,:)/g0*1e6, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, accel_meas_all(1,:)/g0*1e6, 'r-', 'LineWidth', 0.5);
xlabel('Zaman (s)');
ylabel('a_x (µg)');
title('STIM377H Accel - X Ekseni');
legend('Gerçek', 'Ölçüm', 'Location', 'best');
grid on;

% 5. Accelerometer Hataları
subplot(3,4,5);
plot(t, accel_error_mg(1,:), 'r-', 'LineWidth', 0.5);
hold on;
plot(t, accel_error_mg(2,:), 'g-', 'LineWidth', 0.5);
plot(t, accel_error_mg(3,:), 'b-', 'LineWidth', 0.5);
xlabel('Zaman (s)');
ylabel('Hata (mg)');
title('Accelerometer Ölçüm Hataları');
legend('X', 'Y', 'Z', 'Location', 'best');
grid on;

% 6. Magnetometer Ölçümleri
subplot(3,4,6);
plot(t, B_true_all(3,:)/1000, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, B_meas_all(3,:)/1000, 'r-', 'LineWidth', 0.5);
xlabel('Zaman (s)');
ylabel('B_z (µT)');
title('MAG-3 Magnetometer - Z Ekseni');
legend('Gerçek', 'Ölçüm', 'Location', 'best');
grid on;

% 7. Magnetometer Hataları
subplot(3,4,7);
plot(t, mag_error(1,:), 'r-', 'LineWidth', 0.5);
hold on;
plot(t, mag_error(2,:), 'g-', 'LineWidth', 0.5);
plot(t, mag_error(3,:), 'b-', 'LineWidth', 0.5);
xlabel('Zaman (s)');
ylabel('Hata (nT)');
title('Magnetometer Ölçüm Hataları');
legend('X', 'Y', 'Z', 'Location', 'best');
grid on;

% 8. Magnetometer Vektör Büyüklüğü
subplot(3,4,8);
B_mag_true = vecnorm(B_true_all);
B_mag_meas = vecnorm(B_meas_all);
plot(t, B_mag_true/1000, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, B_mag_meas/1000, 'r-', 'LineWidth', 0.5);
xlabel('Zaman (s)');
ylabel('|B| (µT)');
title('Manyetik Alan Büyüklüğü');
legend('Gerçek', 'Ölçüm', 'Location', 'best');
grid on;

% 9. GPS Pozisyon (X-Y düzlemi)
subplot(3,4,9);
plot(r_true_all(1,:), r_true_all(2,:), 'b-', 'LineWidth', 1.5);
hold on;
plot(r_meas_all(1,valid_idx), r_meas_all(2,valid_idx), 'r.', 'MarkerSize', 4);
xlabel('X (km)');
ylabel('Y (km)');
title('GNSS-701 GPS Pozisyon');
legend('Gerçek Yörünge', 'GPS Ölçümler', 'Location', 'best');
axis equal;
grid on;

% 10. GPS Pozisyon Hatası
subplot(3,4,10);
if ~isempty(valid_idx)
    plot(t_gps(valid_idx), vecnorm(pos_error_m), 'b.', 'MarkerSize', 4);
    xlabel('Zaman (s)');
    ylabel('3D Pozisyon Hatası (m)');
    title('GPS Pozisyon Hatası');
    yline(1.5, 'r--', '1.5m spec', 'LineWidth', 1.5);
    grid on;
end

% 11. GPS Hız Hatası
subplot(3,4,11);
if ~isempty(valid_idx)
    plot(t_gps(valid_idx), vecnorm(vel_error_ms), 'b.', 'MarkerSize', 4);
    xlabel('Zaman (s)');
    ylabel('3D Hız Hatası (m/s)');
    title('GPS Hız Hatası');
    yline(0.03, 'r--', '0.03 m/s spec', 'LineWidth', 1.5);
    grid on;
end

% 12. Sensör Özet Tablosu
subplot(3,4,12);
axis off;
text(0.1, 0.95, 'SENSÖR ÖZETİ', 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.80, sprintf('STIM377H Gyro ARW: 0.15 °/√h'), 'FontSize', 9);
text(0.1, 0.70, sprintf('STIM377H Gyro Bias: 0.5 °/h'), 'FontSize', 9);
text(0.1, 0.60, sprintf('STIM377H Accel VRW: 0.07 m/s/√h'), 'FontSize', 9);
text(0.1, 0.50, sprintf('MAG-3 Noise: 12 pT/√Hz'), 'FontSize', 9);
text(0.1, 0.40, sprintf('MAG-3 Accuracy: ±0.75%% FS'), 'FontSize', 9);
text(0.1, 0.30, sprintf('GNSS-701 Pos: 1.5 m'), 'FontSize', 9);
text(0.1, 0.20, sprintf('GNSS-701 Vel: 0.03 m/s'), 'FontSize', 9);

sgtitle('Sensör Modelleri Test Sonuçları', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nTüm testler tamamlandı.\n');
fprintf('========================================\n');
