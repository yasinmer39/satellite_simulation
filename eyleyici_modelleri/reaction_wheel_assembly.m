% 3-Axis Reaction Wheel Assembly Model
% Üç eksenli reaction wheel sistemi

% Bu fonksiyon, 3 eksenli RWA'yı simüle eder.
% Her eksen için ayrı reaction wheel modeli kullanılır.

% inputs:
%   rwa_state   - RWA durumu struct (3 wheel state)
%   torque_cmd  - komut torku body frame'de (N·m) [3x1]
%   dt          - zaman adımı (s)
%   params      - RWA parametreleri (opsiyonel)
%
% outputs:
%   rwa_state   - güncellenmiş RWA durumu
%   torque_actual - uygulanan gerçek tork (N·m) [3x1]
%   h_wheels    - tekerlek açısal momentumu (N·m·s) [3x1]
%   power       - toplam güç tüketimi (W)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: reaction_wheel_model

function [rwa_state, torque_actual, h_wheels, power] = reaction_wheel_assembly(rwa_state, torque_cmd, dt, params)

% Varsayılan parametreler
if nargin < 4 || isempty(params)
    params = rwa_default_params();
end

% RWA durumunu başlat (ilk çağrıda)
if isempty(rwa_state)
    rwa_state = struct();
    rwa_state.wheel = cell(3, 1);
    for i = 1:3
        rwa_state.wheel{i} = [];
    end
end

% Vektörü sütun vektörüne çevir
if isrow(torque_cmd)
    torque_cmd = torque_cmd';
end

% Her tekerlek için torque komutu
% Tekerlekler body frame eksenleri boyunca monte edilmiş varsayılır
torque_actual = zeros(3, 1);
h_wheels = zeros(3, 1);
power = 0;

for i = 1:3
    % Her tekerlek için model çağır
    [rwa_state.wheel{i}, tau, ~] = reaction_wheel_model(...
        rwa_state.wheel{i}, torque_cmd(i), dt, params.wheel_params);
    
    torque_actual(i) = tau;
    h_wheels(i) = rwa_state.wheel{i}.h;
    
    % Güç hesaplama (basitleştirilmiş)
    speed_fraction = abs(rwa_state.wheel{i}.speed) / params.wheel_params.max_speed;
    torque_fraction = abs(torque_cmd(i)) / params.wheel_params.max_torque;
    power = power + params.wheel_params.idle_power + ...
            (torque_fraction + 0.1*speed_fraction) * params.wheel_params.max_power;
end

% Mounting matrix uygulaması (eğer tekerlekler eğik monte edilmişse)
% torque_actual_body = A * torque_actual_wheel
if isfield(params, 'mounting_matrix')
    torque_actual = params.mounting_matrix * torque_actual;
    h_wheels = params.mounting_matrix * h_wheels;
end

end


function params = rwa_default_params()
% Varsayılan RWA parametreleri

% Tek tekerlek parametreleri
params.wheel_params = struct();
params.wheel_params.I_w = 0.005;          % kg·m²
params.wheel_params.synch_speed = 6000;   % rpm
params.wheel_params.max_speed = 6000;     % rpm
params.wheel_params.max_torque = 0.02;    % N·m
params.wheel_params.N0 = 0.025;           % N·m
params.wheel_params.a_param = 0.3;
params.wheel_params.N_c = 1e-4;           % N·m
params.wheel_params.N_c_static = 2e-4;    % N·m
params.wheel_params.f = 1e-7;             % N·m/rpm
params.wheel_params.torque_resolution = 1e-5;  % N·m
params.wheel_params.max_torque_rate = 0.1;     % N·m/s
params.wheel_params.idle_power = 0.5;     % W
params.wheel_params.max_power = 5.0;      % W

% Mounting matrix (3x3 identity for orthogonal mounting)
params.mounting_matrix = eye(3);

% Toplam sistem özellikleri
params.total_momentum_capacity = 3 * params.wheel_params.I_w * ...
    params.wheel_params.max_speed * 2*pi/60;  % N·m·s

end
