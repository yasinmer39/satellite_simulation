% STIM377H IMU - Gyroscope Sensor Model
% Safran STIM377H Datasheet (TS1673 rev.5) parametreleri kullanılmıştır.
% Wertz Chapter 7.8 Gyroscope Models referans alınmıştır.

% Bu fonksiyon, STIM377H gyro ölçümlerini simüle eder.

% omega_true - gerçek açısal hız body frame'de (rad/s) [3x1]
% dt         - zaman adımı (s) - gürültü hesaplaması için
% params     - sensör parametreleri (opsiyonel, varsayılan STIM377H)
% omega_meas - ölçülen açısal hız (rad/s) [3x1]
% gyro_state - güncellenen gyro durumu (bias random walk için)

% STIM377H Gyro Specifications (Table 5-3):
%   Full Scale: ±400 °/s
%   Bias Instability: 0.3 °/h (min), 0.5 °/h (typical) @ 600s
%   Angular Random Walk (ARW): 0.15 °/√h
%   Scale Factor Accuracy: ±500 ppm
%   Non-Linearity: 15 ppm (±200°/s), 20 ppm (±400°/s)
%   Misalignment: 1 mrad
%   Bandwidth: 262 Hz

% Gyro ölçüm modeli (Wertz Eq. 7-134):
%   ω_meas = (1 + k) * ω_true + b + n
%   k: scale factor error
%   b: bias (drift)
%   n: white noise (ARW)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [omega_meas, gyro_state] = gyro_model_stim377h(omega_true, dt, params, gyro_state)

% Varsayılan parametreler (STIM377H datasheet'ten)
if nargin < 3 || isempty(params)
    params = stim377h_gyro_default_params();
end

% Gyro durumunu başlat (ilk çağrıda)
if nargin < 4 || isempty(gyro_state)
    gyro_state = struct();
    gyro_state.bias = params.bias_instability * randn(3,1);  % Başlangıç bias
    gyro_state.bias_rw_sigma = params.bias_instability / sqrt(params.bias_correlation_time);
end

% Vektörü sütun vektörüne çevir
if isrow(omega_true)
    omega_true = omega_true';
end

% 1. Scale Factor Error
% SF accuracy: ±500 ppm = ±0.0005
scale_factor = 1 + params.scale_factor_error;

% 2. Misalignment Error
% Misalignment: 1 mrad
omega_misaligned = params.misalignment_matrix * omega_true;

% 3. Scale factor uygula
omega_scaled = scale_factor .* omega_misaligned;

% 4. Bias (drift) - random walk modeli
% Bias instability update (1st order Gauss-Markov process)
tau = params.bias_correlation_time;
sigma_rw = gyro_state.bias_rw_sigma;

% Gauss-Markov bias dynamics
gyro_state.bias = exp(-dt/tau) * gyro_state.bias + ...
                  sigma_rw * sqrt(1 - exp(-2*dt/tau)) * randn(3,1);

% 5. Angular Random Walk (white noise on rate)
% ARW: 0.15 °/√h = 0.15 * (π/180) / sqrt(3600) rad/√s
arw_sigma = params.arw * sqrt(dt);  % discrete time noise
noise = arw_sigma * randn(3,1);

% 6. Quantization
% Resolution: 0.22 °/h = 1.06e-6 rad/s
omega_quantized = round(omega_scaled / params.resolution) * params.resolution;

% Toplam ölçüm
omega_meas = omega_quantized + gyro_state.bias + noise;

end


function params = stim377h_gyro_default_params()
% STIM377H varsayılan parametreleri

% Birim dönüşümleri
deg2rad = pi/180;
hour2sec = 3600;

% Full scale (saturation için)
params.full_scale = 400 * deg2rad;  % rad/s

% Bias Instability: 0.5 °/h typical
params.bias_instability = 0.5 * deg2rad / hour2sec;  % rad/s

% Angular Random Walk: 0.15 °/√h
params.arw = 0.15 * deg2rad / sqrt(hour2sec);  % rad/√s

% Scale Factor Error: ±500 ppm (rastgele her eksende)
params.scale_factor_error = 500e-6 * (2*rand(3,1) - 1);

% Resolution: 0.22 °/h
params.resolution = 0.22 * deg2rad / hour2sec;  % rad/s

% Misalignment: 1 mrad
misalign_angle = 1e-3;  % rad
% Küçük rastgele misalignment matrisi
theta_x = misalign_angle * (2*rand - 1);
theta_y = misalign_angle * (2*rand - 1);
theta_z = misalign_angle * (2*rand - 1);
params.misalignment_matrix = eye(3) + [0, -theta_z, theta_y;
                                        theta_z, 0, -theta_x;
                                       -theta_y, theta_x, 0];

% Bias correlation time (Gauss-Markov için)
params.bias_correlation_time = 600;  % s (Allan variance integration time)

% Bandwidth
params.bandwidth = 262;  % Hz

% Sample rate
params.sample_rate = 2000;  % Hz

end
