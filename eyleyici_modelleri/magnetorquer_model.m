% Magnetorquer (Magnetic Coil) Actuator Model
% Wertz Chapter 6.7 ve 19.1 referans alınmıştır.
% Denklem 6-20, 19-3

% Bu fonksiyon, magnetorquer davranışını simüle eder.
% Manyetik dipol üretimi ve tork hesaplaması yapılır.

% inputs:
%   dipole_cmd  - komut dipol momenti body frame'de (A·m²) [3x1]
%   B_body      - manyetik alan body frame'de (T) [3x1]
%   params      - magnetorquer parametreleri (opsiyonel)
%
% outputs:
%   dipole_actual - uygulanan gerçek dipol (A·m²) [3x1]
%   torque        - uzay aracına uygulanan tork (N·m) [3x1]
%   power         - güç tüketimi (W)

% Magnetic Dipole (Wertz Eq. 6-20):
%   d = μ * N * I * A * n
%   μ: permeability (air core: μ₀ = 4π×10⁻⁷)
%   N: number of turns
%   I: current
%   A: coil area
%   n: unit normal vector

% Magnetic Torque (Wertz Eq. 19-3):
%   N = M × B = m × B
%   M: magnetic dipole moment
%   B: external magnetic field

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [dipole_actual, torque, power] = magnetorquer_model(dipole_cmd, B_body, params)

% Varsayılan parametreler
if nargin < 3 || isempty(params)
    params = magnetorquer_default_params();
end

% Vektörleri sütun vektörüne çevir
if isrow(dipole_cmd)
    dipole_cmd = dipole_cmd';
end
if isrow(B_body)
    B_body = B_body';
end

% 1. Dipol Saturasyonu
% Her eksen için maksimum dipol kontrolü
dipole_saturated = zeros(3, 1);
for i = 1:3
    dipole_saturated(i) = max(min(dipole_cmd(i), params.max_dipole(i)), -params.max_dipole(i));
end

% 2. Quantization (DAC resolution)
dipole_quantized = round(dipole_saturated / params.dipole_resolution) * params.dipole_resolution;

% 3. Linearity Error
% Gerçek dipol = komut * (1 + linearity_error)
linearity_error = params.linearity * (dipole_quantized / params.max_dipole).^2;
dipole_with_error = dipole_quantized .* (1 + linearity_error);

% 4. Residual Dipole (powered off residual)
dipole_actual = dipole_with_error + params.residual_dipole;

% 5. Cross-axis Coupling
% Küçük çapraz bağlaşım (non-orthogonality)
coupling_matrix = eye(3) + params.coupling_matrix;
dipole_actual = coupling_matrix * dipole_actual;

% 6. Magnetic Torque Calculation (Wertz Eq. 19-3)
% N = M × B
torque = cross(dipole_actual, B_body);

% 7. Power Consumption
% P = Σ (I_i² * R_i) where I = dipole / (N * A)
% Simplified: P ∝ |dipole|² / max_dipole² * max_power
power_fraction = (dipole_quantized ./ params.max_dipole).^2;
power = params.idle_power + sum(power_fraction) * params.max_power_per_axis;

end


function params = magnetorquer_default_params()
% Varsayılan magnetorquer parametreleri
% Tipik küçük uydu magnetorquer değerleri

% Maximum dipole moment per axis (A·m²)
% Typical range: 0.1 - 20 A·m² for small sats
params.max_dipole = [5; 5; 5];  % A·m²

% Dipole resolution (12-bit DAC varsayımı)
params.dipole_resolution = 10 / 4096;  % A·m²

% Linearity error (fraction of command)
params.linearity = 0.01;  % 1%

% Residual magnetic dipole (when off)
params.residual_dipole = [0.001; 0.001; 0.001];  % A·m²

% Cross-axis coupling matrix (small off-diagonal terms)
params.coupling_matrix = [0, 0.01, -0.005;
                          -0.01, 0, 0.008;
                          0.005, -0.008, 0];

% Power consumption
params.idle_power = 0.1;           % W (electronics)
params.max_power_per_axis = 1.0;   % W (at max dipole)

% Coil parameters (for reference)
params.coil_turns = 200;           % Number of turns
params.coil_area = 0.04;           % m² (20cm x 20cm)
params.coil_resistance = 50;       % Ω

% Response time
params.rise_time = 0.01;           % s (electrical time constant)

end
