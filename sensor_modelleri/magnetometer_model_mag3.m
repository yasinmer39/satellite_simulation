% MAG-3 Three Axis Satellite Magnetometer Model
% AAC Clyde Space MAG-3 Datasheet parametreleri kullanılmıştır.
% Wertz Chapter 7.5 Magnetometer Models referans alınmıştır.

% Bu fonksiyon, MAG-3 fluxgate magnetometer ölçümlerini simüle eder.

% B_true     - gerçek manyetik alan body frame'de (nT) [3x1]
% params     - sensör parametreleri (opsiyonel, varsayılan MAG-3)
% mag_state  - sensör durumu (bias drift için, opsiyonel)
% B_meas     - ölçülen manyetik alan (nT) [3x1]

% MAG-3 Specifications (Datasheet):
%   Accuracy: ±0.75% of Full Scale (0.5% typical)
%   Linearity: ±0.015% of Full Scale
%   Sensitivity: 100 µV/nT
%   Range: ±100 µT = ±100,000 nT (±10V output)
%   Noise: 12 pT RMS/√Hz @ 1 Hz
%   Axial Alignment: Better than ±1 degree
%   Zero Shift with Temperature: ±0.6 nT/°C
%   Analog Output @ Zero Field: ±0.025 Volt (±25 nT equivalent)
%   Bandwidth: 500 Hz (-3dB)

% Magnetometer ölçüm modeli (Wertz Eq. 7-79):
%   V = a * (n̂ · H) + V₀
%   veya vektörel: B_meas = A * B_true + b + n

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [B_meas, mag_state] = magnetometer_model_mag3(B_true, params, mag_state)

% Varsayılan parametreler (MAG-3 datasheet'ten)
if nargin < 2 || isempty(params)
    params = mag3_default_params();
end

% Magnetometer durumunu başlat (ilk çağrıda)
if nargin < 3 || isempty(mag_state)
    mag_state = struct();
    mag_state.bias = params.zero_field_offset * (2*rand(3,1) - 1);
    mag_state.temperature = 25;  % °C (referans sıcaklık)
end

% Vektörü sütun vektörüne çevir
if isrow(B_true)
    B_true = B_true';
end

% 1. Misalignment / Calibration Matrix (Wertz Eq. 7-81)
% V = A * H + V₀
% A içerir: scale factor, alignment
B_aligned = params.calibration_matrix * B_true;

% 2. Scale Factor Error (Accuracy: ±0.75% = ±7500 ppm)
scale_factor = 1 + params.scale_factor_error;
B_scaled = scale_factor .* B_aligned;

% 3. Non-linearity
% Linearity: ±0.015% FS = ±15 nT @ ±100,000 nT range
nonlin_error = params.nonlinearity * (B_scaled / params.full_scale).^2 .* sign(B_scaled);
B_nonlin = B_scaled + nonlin_error;

% 4. Zero Field Offset (Bias)
% Analog Output @ Zero Field: ±25 nT equivalent
B_with_bias = B_nonlin + mag_state.bias;

% 5. Temperature Drift
% Zero Shift with Temperature: ±0.6 nT/°C
% (Eğer sıcaklık değişimi simüle edilecekse)
temp_drift = params.temp_coefficient * (mag_state.temperature - 25);
B_with_temp = B_with_bias + temp_drift * ones(3,1);

% 6. Noise
% 12 pT RMS/√Hz @ 1 Hz
% Örnekleme frekansı 100 Hz varsayılan
noise_sigma = params.noise_density * sqrt(params.bandwidth);  % nT RMS
noise = noise_sigma * randn(3,1);

% 7. Quantization (12-bit ADC varsayımı)
% ±100,000 nT / 4096 = ~49 nT resolution
B_quantized = round(B_with_temp / params.resolution) * params.resolution;

% Toplam ölçüm
B_meas = B_quantized + noise;

% Saturation kontrolü
B_meas = max(min(B_meas, params.full_scale), -params.full_scale);

end


function params = mag3_default_params()
% MAG-3 varsayılan parametreleri

% Full scale: ±100 µT = ±100,000 nT
params.full_scale = 100000;  % nT

% Accuracy: ±0.75% FS (scale factor error olarak)
% Her eksen için rastgele
params.scale_factor_error = 0.0075 * (2*rand(3,1) - 1);

% Linearity: ±0.015% FS = ±15 nT
params.nonlinearity = 15;  % nT

% Zero Field Offset: ±25 nT (from ±0.025V @ 100µV/nT)
params.zero_field_offset = 25;  % nT

% Temperature Coefficient: ±0.6 nT/°C
params.temp_coefficient = 0.6;  % nT/°C

% Noise Density: 12 pT/√Hz = 0.012 nT/√Hz
params.noise_density = 0.012;  % nT/√Hz

% Bandwidth: 500 Hz
params.bandwidth = 500;  % Hz

% Sample rate (varsayılan)
params.sample_rate = 100;  % Hz

% Resolution (12-bit ADC varsayımı)
% ±100,000 nT / 2^12 ≈ 48.8 nT
params.resolution = 100000 / 2048;  % nT

% Misalignment / Orthogonality: Better than ±1 degree
misalign_angle = 1 * pi/180;  % 1 degree in radians
theta_x = misalign_angle * (2*rand - 1);
theta_y = misalign_angle * (2*rand - 1);
theta_z = misalign_angle * (2*rand - 1);

% Calibration matrix (Wertz Eq. 7-81: V = A*H + V0)
% Küçük açılar için: A ≈ I + skew(θ)
params.calibration_matrix = eye(3) + [0, -theta_z, theta_y;
                                       theta_z, 0, -theta_x;
                                      -theta_y, theta_x, 0];

% Perming susceptibility: ±8 nT shift with ±5 Gauss applied
params.perming_susceptibility = 8;  % nT

end
