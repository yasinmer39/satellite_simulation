% STIM377H IMU - Accelerometer Sensor Model
% Safran STIM377H Datasheet (TS1673 rev.5) parametreleri kullanılmıştır.
% 10g accelerometer konfigürasyonu (Table 5-5)

% Bu fonksiyon, STIM377H accelerometer ölçümlerini simüle eder.

% accel_true - gerçek ivme body frame'de (m/s²) [3x1]
% dt         - zaman adımı (s) - gürültü hesaplaması için
% params     - sensör parametreleri (opsiyonel, varsayılan STIM377H 10g)
% accel_meas - ölçülen ivme (m/s²) [3x1]
% accel_state - güncellenen accelerometer durumu (bias random walk için)

% STIM377H 10g Accelerometer Specifications (Table 5-5):
%   Full Scale: ±10 g
%   Resolution: 1.9 µg
%   Bias Instability: 0.04 mg (min), 0.05 mg (typical) @ 600s
%   Velocity Random Walk (VRW): 0.07 m/s/√h
%   Scale Factor Accuracy: ±200 ppm
%   Non-Linearity: 100 ppm
%   Misalignment: 1 mrad
%   Bandwidth: 214 Hz (typical)

% Accelerometer ölçüm modeli:
%   a_meas = (1 + k) * a_true + b + n
%   k: scale factor error
%   b: bias
%   n: white noise (VRW)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [accel_meas, accel_state] = accel_model_stim377h(accel_true, dt, params, accel_state)

% Sabitler
g0 = 9.80665;  % m/s² (standard gravity)

% Varsayılan parametreler (STIM377H 10g datasheet'ten)
if nargin < 3 || isempty(params)
    params = stim377h_accel_default_params();
end

% Accelerometer durumunu başlat (ilk çağrıda)
if nargin < 4 || isempty(accel_state)
    accel_state = struct();
    accel_state.bias = params.bias_instability * randn(3,1);  % Başlangıç bias
    accel_state.bias_rw_sigma = params.bias_instability / sqrt(params.bias_correlation_time);
end

% Vektörü sütun vektörüne çevir
if isrow(accel_true)
    accel_true = accel_true';
end

% 1. Scale Factor Error
% SF accuracy: ±200 ppm
scale_factor = 1 + params.scale_factor_error;

% 2. Misalignment Error
% Misalignment: 1 mrad
accel_misaligned = params.misalignment_matrix * accel_true;

% 3. Scale factor uygula
accel_scaled = scale_factor .* accel_misaligned;

% 4. Bias - random walk modeli
% Bias instability update (1st order Gauss-Markov process)
tau = params.bias_correlation_time;
sigma_rw = accel_state.bias_rw_sigma;

% Gauss-Markov bias dynamics
accel_state.bias = exp(-dt/tau) * accel_state.bias + ...
                   sigma_rw * sqrt(1 - exp(-2*dt/tau)) * randn(3,1);

% 5. Velocity Random Walk (white noise on acceleration)
% VRW: 0.07 m/s/√h
vrw_sigma = params.vrw * sqrt(dt);  % discrete time noise
noise = vrw_sigma * randn(3,1);

% 6. Quantization
% Resolution: 1.9 µg = 1.9e-6 * 9.80665 m/s²
accel_quantized = round(accel_scaled / params.resolution) * params.resolution;

% Toplam ölçüm
accel_meas = accel_quantized + accel_state.bias + noise;

% Saturation kontrolü
accel_meas = max(min(accel_meas, params.full_scale), -params.full_scale);

end


function params = stim377h_accel_default_params()
% STIM377H 10g accelerometer varsayılan parametreleri

% Sabitler
g0 = 9.80665;  % m/s²
hour2sec = 3600;

% Full scale: ±10 g
params.full_scale = 10 * g0;  % m/s²

% Bias Instability: 0.05 mg typical
params.bias_instability = 0.05e-3 * g0;  % m/s²

% Velocity Random Walk: 0.07 m/s/√h
params.vrw = 0.07 / sqrt(hour2sec);  % m/s/√s

% Scale Factor Error: ±200 ppm (rastgele her eksende)
params.scale_factor_error = 200e-6 * (2*rand(3,1) - 1);

% Resolution: 1.9 µg
params.resolution = 1.9e-6 * g0;  % m/s²

% Misalignment: 1 mrad
misalign_angle = 1e-3;  % rad
theta_x = misalign_angle * (2*rand - 1);
theta_y = misalign_angle * (2*rand - 1);
theta_z = misalign_angle * (2*rand - 1);
params.misalignment_matrix = eye(3) + [0, -theta_z, theta_y;
                                        theta_z, 0, -theta_x;
                                       -theta_y, theta_x, 0];

% Bias correlation time (Gauss-Markov için)
params.bias_correlation_time = 600;  % s

% Bandwidth
params.bandwidth = 214;  % Hz

% Sample rate
params.sample_rate = 2000;  % Hz

% Bias error over temperature (for reference)
params.bias_temp_coeff = 1.5e-3 * g0;  % m/s² rms over temp range

end
