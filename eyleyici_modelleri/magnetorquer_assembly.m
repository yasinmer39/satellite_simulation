% 3-Axis Magnetorquer Assembly Model
% Üç eksenli magnetorquer sistemi

% Bu fonksiyon, 3 eksenli magnetorquer sistemini simüle eder.
% B-dot kontrol algoritması opsiyonel olarak dahildir.

% inputs:
%   dipole_cmd  - komut dipol momenti body frame'de (A·m²) [3x1]
%                 veya 'bdot' string'i ile B_body verilirse otomatik B-dot
%   B_body      - manyetik alan body frame'de (T) [3x1]
%   params      - magnetorquer parametreleri (opsiyonel)
%   control_mode - 'direct' veya 'bdot' (opsiyonel)
%   B_body_prev - önceki manyetik alan (B-dot için) [3x1]
%   dt          - zaman adımı (B-dot için) (s)
%
% outputs:
%   dipole_actual - uygulanan gerçek dipol (A·m²) [3x1]
%   torque        - uzay aracına uygulanan tork (N·m) [3x1]
%   power         - güç tüketimi (W)

% B-dot Control Law:
%   M = -k * dB/dt
%   Simple detumbling algorithm

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: magnetorquer_model

function [dipole_actual, torque, power] = magnetorquer_assembly(dipole_cmd, B_body, params, control_mode, B_body_prev, dt)

% Varsayılan parametreler
if nargin < 3 || isempty(params)
    params = mtq_assembly_default_params();
end

if nargin < 4 || isempty(control_mode)
    control_mode = 'direct';
end

% Vektörleri sütun vektörüne çevir
if isrow(B_body)
    B_body = B_body';
end

% B-dot control mode
if strcmp(control_mode, 'bdot')
    if nargin < 5 || isempty(B_body_prev)
        error('B_body_prev required for bdot mode');
    end
    if nargin < 6 || isempty(dt)
        dt = 0.1;  % Default 10 Hz
    end
    
    if isrow(B_body_prev)
        B_body_prev = B_body_prev';
    end
    
    % B-dot hesaplama
    B_dot = (B_body - B_body_prev) / dt;
    
    % B-dot control law: M = -k * dB/dt
    dipole_cmd = -params.bdot_gain * B_dot;
end

% Vektörü sütun vektörüne çevir
if isrow(dipole_cmd)
    dipole_cmd = dipole_cmd';
end

% Magnetorquer model çağır
[dipole_actual, torque, power] = magnetorquer_model(dipole_cmd, B_body, params.mtq_params);

end


function params = mtq_assembly_default_params()
% Varsayılan magnetorquer assembly parametreleri

% Tek magnetorquer parametreleri
params.mtq_params = struct();
params.mtq_params.max_dipole = [5; 5; 5];  % A·m²
params.mtq_params.dipole_resolution = 10 / 4096;  % A·m²
params.mtq_params.linearity = 0.01;
params.mtq_params.residual_dipole = [0.001; 0.001; 0.001];  % A·m²
params.mtq_params.coupling_matrix = [0, 0.01, -0.005;
                                     -0.01, 0, 0.008;
                                     0.005, -0.008, 0];
params.mtq_params.idle_power = 0.1;        % W
params.mtq_params.max_power_per_axis = 1.0; % W

% B-dot control gain
% k should be chosen based on:
%   - Spacecraft inertia
%   - Desired detumbling time
%   - Available dipole moment
% Typical: k = 2 * (1 + sin(inclination)) * sqrt(I / B_max)
params.bdot_gain = 1e6;  % A·m²·s/T

end
