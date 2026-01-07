% Statistical Comparison: True vs Estimated Outputs
% Case Study (e): Gerçek çıktılar ile kestirim çıktılarını istatistiksel olarak karşılaştırınız
%
% Bu script, attitude estimation algoritmalarının performansını
% kapsamlı istatistiksel metriklerle analiz eder.
%
% İstatistiksel Metrikler:
%   - Mean Error (Bias)
%   - Standard Deviation (Precision)
%   - Root Mean Square Error (RMSE)
%   - Maximum/Minimum Error
%   - 3-Sigma Bounds (99.7% confidence)
%   - Histogram ve PDF analizi
%   - Autocorrelation (hata bağımsızlığı)
%   - Normality test (Gaussian assumption)
%
% Referanslar:
%   - Wertz Chapter 12.3 "Covariance Analysis"
%   - Crassidis & Junkins "Optimal Estimation of Dynamic Systems"

clear; clc; close all;

fprintf('=================================================================\n');
fprintf('   CASE STUDY (e): İSTATİSTİKSEL KARŞILAŞTIRMA ANALİZİ\n');
fprintf('   Gerçek Çıktılar vs Kestirim Çıktıları\n');
fprintf('=================================================================\n\n');

%% ==================== SIMULATION SETUP ====================
fprintf('1. Simülasyon Hazırlanıyor...\n');

% Simülasyon parametreleri
t_sim = 1800;        % 30 dakika (daha uzun simülasyon)
dt = 0.1;            % 10 Hz
t = 0:dt:t_sim;
N = length(t);

% Başlangıç zamanı
utc_start = datetime(2026, 6, 21, 12, 0, 0);

% Yörünge parametreleri
mu = 3.986004418e14;
R_earth = 6378137;
altitude = 1000e3;
a = R_earth + altitude;
n_orbit = sqrt(mu / a^3);
inc = 60 * pi/180;

%% ==================== TRUE STATE GENERATION ====================
fprintf('2. Gerçek Durum (Ground Truth) Oluşturuluyor...\n');

% Rastgele başlangıç yönelimi
rng(42);  % Reproducibility için
euler_init = [15; -8; 45] * pi/180;

% Değişken angular velocity (daha realistik)
omega_base = [0.02; -0.01; 0.015] * pi/180;  % rad/s

% True quaternion trajectory
q_true = zeros(4, N);
q_true(:,1) = euler_to_quat_stat(euler_init);

omega_true = zeros(3, N);
for k = 1:N
    % Yavaş değişen angular velocity
    omega_true(:,k) = omega_base .* (1 + 0.1*sin(2*pi*t(k)/600));
end

for k = 2:N
    q_true(:,k) = propagate_quat_stat(q_true(:,k-1), omega_true(:,k-1), dt);
end

%% ==================== ENVIRONMENT MODELS ====================
fprintf('3. Çevre Modelleri Hesaplanıyor...\n');

% Position, magnetic field, sun, gravity
r_eci = zeros(3, N);
B_eci = zeros(3, N);
s_eci = zeros(3, N);
g_eci = zeros(3, N);

for k = 1:N
    % Orbital position
    M = n_orbit * t(k);
    r_orbital = a * [cos(M); sin(M); 0];
    R_inc = [1, 0, 0; 0, cos(inc), -sin(inc); 0, sin(inc), cos(inc)];
    r_eci(:,k) = R_inc * r_orbital;
    
    % Sun position
    utc_now = utc_start + seconds(t(k));
    s_eci(:,k) = sun_pos_stat(utc_now);
    
    % Magnetic field (dipole)
    B_eci(:,k) = mag_field_stat(r_eci(:,k));
    
    % Gravity direction
    g_eci(:,k) = -r_eci(:,k) / norm(r_eci(:,k));
end

%% ==================== SENSOR MEASUREMENTS ====================
fprintf('4. Sensör Ölçümleri Simüle Ediliyor...\n');

% Sensor noise parameters (STIM377H, MAG-3 datasheet'lerden)
gyro_arw = 0.15 * (pi/180) / sqrt(3600);      % 0.15 deg/sqrt(h)
gyro_bias_inst = 0.5 * (pi/180) / 3600;       % 0.5 deg/h
gyro_bias_true = [0.4; -0.3; 0.5] * (pi/180) / 3600;

accel_noise = 100e-6 * 9.81;                   % 100 μg
mag_noise_pct = 0.02;                          % 2% of measurement

% Generate measurements
gyro_meas = zeros(3, N);
accel_meas = zeros(3, N);
mag_meas = zeros(3, N);

for k = 1:N
    A_true = quat_to_dcm_stat(q_true(:,k));
    
    % Gyro: true rate + bias + ARW noise
    gyro_meas(:,k) = omega_true(:,k) + gyro_bias_true + ...
                     gyro_arw * sqrt(1/dt) * randn(3,1);
    
    % Accelerometer: gravity in body + noise
    g_body = A_true * g_eci(:,k);
    accel_meas(:,k) = -g_body * 9.81 + accel_noise * randn(3,1);
    
    % Magnetometer: B field in body + noise
    B_body = A_true * B_eci(:,k);
    mag_meas(:,k) = B_body .* (1 + mag_noise_pct * randn(3,1));
end

%% ==================== ESTIMATION ALGORITHMS ====================
fprintf('5. Kestirim Algoritmaları Çalıştırılıyor...\n');

% Storage for estimates
q_triad = zeros(4, N);
q_quest = zeros(4, N);
q_mekf = zeros(4, N);
bias_mekf = zeros(3, N);

% MEKF initialization
mekf_state = [];
P_mekf = [];
mekf_params = struct();
mekf_params.gyro_noise = gyro_arw;
mekf_params.gyro_bias_noise = gyro_bias_inst;
mekf_params.accel_noise = 0.02;
mekf_params.mag_noise = 0.03;
mekf_params.init_attitude_var = (15*pi/180)^2;
mekf_params.init_bias_var = (1*pi/180/3600)^2;

for k = 1:N
    % Normalize measurements
    accel_norm = -accel_meas(:,k) / norm(accel_meas(:,k));
    mag_norm = mag_meas(:,k) / norm(mag_meas(:,k));
    
    % Reference vectors
    g_ref = g_eci(:,k);
    B_ref = B_eci(:,k) / norm(B_eci(:,k));
    
    % TRIAD
    [q_triad(:,k), ~] = triad_algorithm(accel_norm, mag_norm, g_ref, B_ref);
    
    % QUEST
    v_body = [accel_norm, mag_norm];
    v_ref = [g_ref, B_ref];
    [q_quest(:,k), ~, ~] = quest_algorithm(v_body, v_ref, [0.7, 0.3]);
    
    % MEKF
    [q_mekf(:,k), mekf_state, P_mekf] = mekf_attitude(...
        gyro_meas(:,k), accel_norm, mag_norm, B_ref, g_ref, dt, ...
        mekf_state, P_mekf, mekf_params);
    bias_mekf(:,k) = mekf_state.gyro_bias;
end

%% ==================== ERROR CALCULATION ====================
fprintf('6. Hata Hesaplamaları...\n');

% Attitude errors (degrees)
err_triad = zeros(1, N);
err_quest = zeros(1, N);
err_mekf = zeros(1, N);

% Per-axis errors (Euler angles)
euler_err_triad = zeros(3, N);
euler_err_quest = zeros(3, N);
euler_err_mekf = zeros(3, N);

% Sun vector errors
sun_err_mekf = zeros(1, N);

for k = 1:N
    % Total attitude error
    err_triad(k) = quat_angle_err_stat(q_triad(:,k), q_true(:,k));
    err_quest(k) = quat_angle_err_stat(q_quest(:,k), q_true(:,k));
    err_mekf(k) = quat_angle_err_stat(q_mekf(:,k), q_true(:,k));
    
    % Per-axis Euler errors
    euler_true = quat_to_euler_stat(q_true(:,k)) * 180/pi;
    euler_triad = quat_to_euler_stat(q_triad(:,k)) * 180/pi;
    euler_quest = quat_to_euler_stat(q_quest(:,k)) * 180/pi;
    euler_mekf = quat_to_euler_stat(q_mekf(:,k)) * 180/pi;
    
    euler_err_triad(:,k) = euler_triad - euler_true;
    euler_err_quest(:,k) = euler_quest - euler_true;
    euler_err_mekf(:,k) = euler_mekf - euler_true;
    
    % Sun vector error
    A_true = quat_to_dcm_stat(q_true(:,k));
    A_est = quat_to_dcm_stat(q_mekf(:,k));
    s_body_true = A_true * s_eci(:,k);
    s_body_est = A_est * s_eci(:,k);
    cos_sun = dot(s_body_true, s_body_est);
    sun_err_mekf(k) = acosd(min(1, max(-1, cos_sun)));
end

% Gyro bias error
bias_err = (bias_mekf - gyro_bias_true) * 3600 * 180/pi;  % deg/h

%% ==================== STATISTICAL METRICS ====================
fprintf('\n');
fprintf('=================================================================\n');
fprintf('                    İSTATİSTİKSEL ANALİZ SONUÇLARI\n');
fprintf('=================================================================\n\n');

% Function to compute statistics
compute_stats = @(x) struct('mean', mean(x), 'std', std(x), ...
    'rms', sqrt(mean(x.^2)), 'max', max(x), 'min', min(x), ...
    'p95', prctile(x, 95), 'p99', prctile(x, 99), ...
    'sigma3', mean(x) + 3*std(x));

% Compute statistics for each algorithm
stats_triad = compute_stats(err_triad);
stats_quest = compute_stats(err_quest);
stats_mekf = compute_stats(err_mekf);
stats_sun = compute_stats(sun_err_mekf);

fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('                    ATTITUDE ERROR COMPARISON\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('%-15s %12s %12s %12s\n', 'Metric', 'TRIAD', 'QUEST', 'MEKF');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', 'Mean (Bias)', stats_triad.mean, stats_quest.mean, stats_mekf.mean);
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', 'Std Dev', stats_triad.std, stats_quest.std, stats_mekf.std);
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', 'RMSE', stats_triad.rms, stats_quest.rms, stats_mekf.rms);
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', 'Maximum', stats_triad.max, stats_quest.max, stats_mekf.max);
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', 'Minimum', stats_triad.min, stats_quest.min, stats_mekf.min);
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', '95th Pctl', stats_triad.p95, stats_quest.p95, stats_mekf.p95);
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', '99th Pctl', stats_triad.p99, stats_quest.p99, stats_mekf.p99);
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', '3σ Bound', stats_triad.sigma3, stats_quest.sigma3, stats_mekf.sigma3);
fprintf('─────────────────────────────────────────────────────────────────\n');

fprintf('\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('                    PER-AXIS EULER ANGLE ERRORS (MEKF)\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('%-15s %12s %12s %12s\n', 'Metric', 'Roll', 'Pitch', 'Yaw');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', 'Mean', mean(euler_err_mekf(1,:)), mean(euler_err_mekf(2,:)), mean(euler_err_mekf(3,:)));
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', 'Std Dev', std(euler_err_mekf(1,:)), std(euler_err_mekf(2,:)), std(euler_err_mekf(3,:)));
fprintf('%-15s %12.4f %12.4f %12.4f deg\n', 'RMSE', rms(euler_err_mekf(1,:)), rms(euler_err_mekf(2,:)), rms(euler_err_mekf(3,:)));
fprintf('─────────────────────────────────────────────────────────────────\n');

fprintf('\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('                    GYRO BIAS ESTIMATION ERROR\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('%-15s %12s %12s %12s\n', 'Metric', 'bx', 'by', 'bz');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('%-15s %12.4f %12.4f %12.4f deg/h\n', 'True Bias', gyro_bias_true(1)*3600*180/pi, gyro_bias_true(2)*3600*180/pi, gyro_bias_true(3)*3600*180/pi);
fprintf('%-15s %12.4f %12.4f %12.4f deg/h\n', 'Final Est', bias_mekf(1,end)*3600*180/pi, bias_mekf(2,end)*3600*180/pi, bias_mekf(3,end)*3600*180/pi);
fprintf('%-15s %12.4f %12.4f %12.4f deg/h\n', 'Final Error', bias_err(1,end), bias_err(2,end), bias_err(3,end));
fprintf('%-15s %12.4f deg/h\n', 'Total Error', norm(bias_err(:,end)));
fprintf('─────────────────────────────────────────────────────────────────\n');

fprintf('\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('                    SUN VECTOR ESTIMATION ERROR\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('%-15s %12.4f deg\n', 'Mean', stats_sun.mean);
fprintf('%-15s %12.4f deg\n', 'Std Dev', stats_sun.std);
fprintf('%-15s %12.4f deg\n', 'RMSE', stats_sun.rms);
fprintf('%-15s %12.4f deg\n', 'Maximum', stats_sun.max);
fprintf('%-15s %12.4f deg\n', '3σ Bound', stats_sun.sigma3);
fprintf('─────────────────────────────────────────────────────────────────\n');

%% ==================== NORMALITY TEST ====================
fprintf('\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('                    NORMALITY TEST (Lilliefors)\n');
fprintf('─────────────────────────────────────────────────────────────────\n');

% Test if errors are normally distributed
try
    [h_triad, p_triad] = lillietest(err_triad);
    [h_quest, p_quest] = lillietest(err_quest);
    [h_mekf, p_mekf] = lillietest(err_mekf);
    
    fprintf('%-15s %12s %12s\n', 'Algorithm', 'p-value', 'Normal?');
    fprintf('─────────────────────────────────────────────────────────────────\n');
    fprintf('%-15s %12.4f %12s\n', 'TRIAD', p_triad, ternary(h_triad==0, 'YES', 'NO'));
    fprintf('%-15s %12.4f %12s\n', 'QUEST', p_quest, ternary(h_quest==0, 'YES', 'NO'));
    fprintf('%-15s %12.4f %12s\n', 'MEKF', p_mekf, ternary(h_mekf==0, 'YES', 'NO'));
catch
    fprintf('   Lillietest not available. Skipping normality test.\n');
end
fprintf('─────────────────────────────────────────────────────────────────\n');

%% ==================== AUTOCORRELATION ====================
fprintf('\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('                    AUTOCORRELATION (lag-1)\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('   (Ideal: ~0 for independent errors)\n\n');

acf_triad = autocorr_stat(err_triad, 1);
acf_quest = autocorr_stat(err_quest, 1);
acf_mekf = autocorr_stat(err_mekf, 1);

fprintf('%-15s %12.4f\n', 'TRIAD', acf_triad);
fprintf('%-15s %12.4f\n', 'QUEST', acf_quest);
fprintf('%-15s %12.4f\n', 'MEKF', acf_mekf);
fprintf('─────────────────────────────────────────────────────────────────\n');

%% ==================== CONVERGENCE ANALYSIS ====================
fprintf('\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('                    MEKF CONVERGENCE ANALYSIS\n');
fprintf('─────────────────────────────────────────────────────────────────\n');

% Find convergence time (when error < 2*final_std)
final_std = std(err_mekf(end-1000:end));
conv_threshold = 2 * final_std + mean(err_mekf(end-1000:end));
conv_idx = find(err_mekf < conv_threshold, 1, 'first');
conv_time = t(conv_idx);

fprintf('%-25s %10.1f s\n', 'Convergence Time', conv_time);
fprintf('%-25s %10.4f deg\n', 'Steady-State Mean', mean(err_mekf(end-1000:end)));
fprintf('%-25s %10.4f deg\n', 'Steady-State Std', final_std);
fprintf('%-25s %10.4f deg\n', 'Initial Error', err_mekf(1));
fprintf('─────────────────────────────────────────────────────────────────\n');

%% ==================== PLOTS ====================
fprintf('\n7. Grafikler Oluşturuluyor...\n');

% Figure 1: Error Time Series
figure('Name', 'Attitude Error Comparison', 'Position', [50, 50, 1400, 900]);

subplot(3,2,1);
plot(t, err_triad, 'b-', 'LineWidth', 0.5);
hold on;
plot(t, err_quest, 'r-', 'LineWidth', 0.5);
plot(t, err_mekf, 'g-', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Error (deg)');
title('Attitude Estimation Error vs Time');
legend('TRIAD', 'QUEST', 'MEKF', 'Location', 'best');
grid on;

subplot(3,2,2);
plot(t, err_mekf, 'g-', 'LineWidth', 1);
hold on;
yline(stats_mekf.mean, 'k--', 'Mean');
yline(stats_mekf.sigma3, 'r--', '3\sigma');
xlabel('Time (s)'); ylabel('MEKF Error (deg)');
title('MEKF Error with Statistical Bounds');
grid on;
legend('Error', 'Mean', '3\sigma Bound');

subplot(3,2,3);
histogram(err_triad, 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5);
hold on;
histogram(err_quest, 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5);
histogram(err_mekf, 50, 'Normalization', 'pdf', 'FaceAlpha', 0.5);
xlabel('Error (deg)'); ylabel('Probability Density');
title('Error Distribution (PDF)');
legend('TRIAD', 'QUEST', 'MEKF');
grid on;

subplot(3,2,4);
plot(t, euler_err_mekf(1,:), 'r-', 'LineWidth', 0.8);
hold on;
plot(t, euler_err_mekf(2,:), 'g-', 'LineWidth', 0.8);
plot(t, euler_err_mekf(3,:), 'b-', 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('Error (deg)');
title('MEKF Per-Axis Euler Angle Error');
legend('Roll', 'Pitch', 'Yaw');
grid on;

subplot(3,2,5);
plot(t, bias_mekf(1,:)*3600*180/pi, 'r-', 'LineWidth', 1);
hold on;
plot(t, bias_mekf(2,:)*3600*180/pi, 'g-', 'LineWidth', 1);
plot(t, bias_mekf(3,:)*3600*180/pi, 'b-', 'LineWidth', 1);
yline(gyro_bias_true(1)*3600*180/pi, 'r--');
yline(gyro_bias_true(2)*3600*180/pi, 'g--');
yline(gyro_bias_true(3)*3600*180/pi, 'b--');
xlabel('Time (s)'); ylabel('Bias (deg/h)');
title('Gyro Bias Estimation');
legend('b_x est', 'b_y est', 'b_z est', 'b_x true', 'b_y true', 'b_z true');
grid on;

subplot(3,2,6);
plot(t, sun_err_mekf, 'm-', 'LineWidth', 1);
hold on;
yline(stats_sun.mean, 'k--');
yline(stats_sun.sigma3, 'r--');
xlabel('Time (s)'); ylabel('Error (deg)');
title('Sun Vector Estimation Error (Body Frame)');
legend('Error', 'Mean', '3\sigma');
grid on;

saveas(gcf, 'statistical_comparison_timeseries.png');

% Figure 2: Box plots and QQ plots
figure('Name', 'Statistical Analysis', 'Position', [100, 100, 1200, 600]);

subplot(2,3,1);
boxplot([err_triad', err_quest', err_mekf'], {'TRIAD', 'QUEST', 'MEKF'});
ylabel('Attitude Error (deg)');
title('Error Distribution (Box Plot)');
grid on;

subplot(2,3,2);
qqplot(err_mekf);
title('MEKF Error Q-Q Plot');
grid on;

subplot(2,3,3);
% Autocorrelation plot
[acf_vals, lags] = autocorr_stat_full(err_mekf, 50);
stem(lags, acf_vals, 'filled');
xlabel('Lag'); ylabel('Autocorrelation');
title('MEKF Error Autocorrelation');
grid on;

subplot(2,3,4);
% Cumulative distribution
[f_triad, x_triad] = ecdf(err_triad);
[f_quest, x_quest] = ecdf(err_quest);
[f_mekf, x_mekf] = ecdf(err_mekf);
plot(x_triad, f_triad, 'b-', 'LineWidth', 1.5);
hold on;
plot(x_quest, f_quest, 'r-', 'LineWidth', 1.5);
plot(x_mekf, f_mekf, 'g-', 'LineWidth', 1.5);
xlabel('Error (deg)'); ylabel('CDF');
title('Cumulative Distribution Function');
legend('TRIAD', 'QUEST', 'MEKF', 'Location', 'southeast');
grid on;

subplot(2,3,5);
% Error vs time (moving average)
window = 100;
ma_triad = movmean(err_triad, window);
ma_quest = movmean(err_quest, window);
ma_mekf = movmean(err_mekf, window);
plot(t, ma_triad, 'b-', 'LineWidth', 1);
hold on;
plot(t, ma_quest, 'r-', 'LineWidth', 1);
plot(t, ma_mekf, 'g-', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Error (deg)');
title(['Moving Average (window=' num2str(window) ')']);
legend('TRIAD', 'QUEST', 'MEKF');
grid on;

subplot(2,3,6);
% Bias estimation error
plot(t, vecnorm(bias_err), 'k-', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Bias Error (deg/h)');
title('Total Gyro Bias Estimation Error');
grid on;

saveas(gcf, 'statistical_comparison_analysis.png');

fprintf('   Grafikler kaydedildi:\n');
fprintf('   - statistical_comparison_timeseries.png\n');
fprintf('   - statistical_comparison_analysis.png\n');

%% ==================== SUMMARY TABLE (LaTeX format) ====================
fprintf('\n');
fprintf('=================================================================\n');
fprintf('                    ÖZET TABLO (Rapor için)\n');
fprintf('=================================================================\n\n');

fprintf('\\begin{table}[h]\n');
fprintf('\\centering\n');
fprintf('\\caption{Attitude Estimation Algorithm Comparison}\n');
fprintf('\\begin{tabular}{|l|c|c|c|}\n');
fprintf('\\hline\n');
fprintf('Metric & TRIAD & QUEST & MEKF \\\\ \\hline\n');
fprintf('Mean Error (deg) & %.4f & %.4f & %.4f \\\\ \\hline\n', stats_triad.mean, stats_quest.mean, stats_mekf.mean);
fprintf('Std Dev (deg) & %.4f & %.4f & %.4f \\\\ \\hline\n', stats_triad.std, stats_quest.std, stats_mekf.std);
fprintf('RMSE (deg) & %.4f & %.4f & %.4f \\\\ \\hline\n', stats_triad.rms, stats_quest.rms, stats_mekf.rms);
fprintf('3$\\sigma$ Bound (deg) & %.4f & %.4f & %.4f \\\\ \\hline\n', stats_triad.sigma3, stats_quest.sigma3, stats_mekf.sigma3);
fprintf('\\end{tabular}\n');
fprintf('\\end{table}\n');

fprintf('\n=== ANALİZ TAMAMLANDI ===\n');

%% ==================== LOCAL HELPER FUNCTIONS ====================

function q = euler_to_quat_stat(e)
    cr=cos(e(1)/2); sr=sin(e(1)/2);
    cp=cos(e(2)/2); sp=sin(e(2)/2);
    cy=cos(e(3)/2); sy=sin(e(3)/2);
    q = [sr*cp*cy-cr*sp*sy; cr*sp*cy+sr*cp*sy; cr*cp*sy-sr*sp*cy; cr*cp*cy+sr*sp*sy];
    q = q/norm(q);
end

function e = quat_to_euler_stat(q)
    q1=q(1); q2=q(2); q3=q(3); q4=q(4);
    e = [atan2(2*(q4*q1+q2*q3), 1-2*(q1^2+q2^2));
         asin(max(-1,min(1,2*(q4*q2-q3*q1))));
         atan2(2*(q4*q3+q1*q2), 1-2*(q2^2+q3^2))];
end

function q = propagate_quat_stat(q, w, dt)
    Omega = [0,w(3),-w(2),w(1); -w(3),0,w(1),w(2); w(2),-w(1),0,w(3); -w(1),-w(2),-w(3),0];
    q = q + 0.5*Omega*q*dt;
    q = q/norm(q);
end

function A = quat_to_dcm_stat(q)
    q1=q(1); q2=q(2); q3=q(3); q4=q(4);
    A = [q1^2-q2^2-q3^2+q4^2, 2*(q1*q2+q3*q4), 2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4), -q1^2+q2^2-q3^2+q4^2, 2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4), 2*(q2*q3-q1*q4), -q1^2-q2^2+q3^2+q4^2];
end

function err = quat_angle_err_stat(q1, q2)
    dq = quat_mult_stat(q2, [-q1(1:3); q1(4)]);
    if dq(4)<0, dq=-dq; end
    err = 2*acos(min(1,abs(dq(4))))*180/pi;
end

function r = quat_mult_stat(p, q)
    r = [p(4)*q(1)+p(1)*q(4)+p(2)*q(3)-p(3)*q(2);
         p(4)*q(2)-p(1)*q(3)+p(2)*q(4)+p(3)*q(1);
         p(4)*q(3)+p(1)*q(2)-p(2)*q(1)+p(3)*q(4);
         p(4)*q(4)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)];
end

function s = sun_pos_stat(utc)
    JD = juliandate(utc);
    T = (JD-2451545)/36525;
    L = mod(280.46+36000.77*T, 360);
    M = mod(357.53+35999.05*T, 360)*pi/180;
    lam = (L + 1.915*sin(M) + 0.02*sin(2*M))*pi/180;
    eps = (23.439-0.013*T)*pi/180;
    s = [cos(lam); cos(eps)*sin(lam); sin(eps)*sin(lam)];
    s = s/norm(s);
end

function B = mag_field_stat(r)
    B0 = 3.12e-5; Re = 6378137;
    r_mag = norm(r);
    r_hat = r/r_mag;
    m_hat = [0;0;1];
    B = B0*(Re/r_mag)^3*(3*dot(m_hat,r_hat)*r_hat - m_hat);
end

function acf = autocorr_stat(x, lag)
    x = x - mean(x);
    n = length(x);
    acf = sum(x(1:n-lag).*x(1+lag:n)) / sum(x.^2);
end

function [acf, lags] = autocorr_stat_full(x, maxlag)
    lags = 0:maxlag;
    acf = zeros(1, maxlag+1);
    for k = 0:maxlag
        acf(k+1) = autocorr_stat(x, k);
    end
end

function r = ternary(cond, a, b)
    if cond, r = a; else, r = b; end
end
