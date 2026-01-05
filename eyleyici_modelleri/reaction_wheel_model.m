% Reaction Wheel Actuator Model
% Wertz Chapter 7.9 - Reaction Wheel Models referans alınmıştır.
% Denklem 7-144, 7-145, 7-146

% Bu fonksiyon, reaction wheel davranışını simüle eder.
% Tekerlek torku, hız ve sürtünme modellenir.

% inputs:
%   wheel_state - tekerlek durumu struct:
%                 .speed - tekerlek hızı (rpm)
%                 .h     - açısal momentum (N·m·s)
%   torque_cmd  - komut torku (N·m) [skaler]
%   dt          - zaman adımı (s)
%   params      - tekerlek parametreleri (opsiyonel)
%
% outputs:
%   wheel_state - güncellenmiş tekerlek durumu
%   torque_actual - uygulanan gerçek tork (N·m)
%   h_dot       - açısal momentum türevi (N·m)

% Reaction Wheel Dynamics (Wertz Eq. 7-144):
%   N = X_dc * N_em - N_friction
%   X_dc: duty cycle (-1 to +1)
%   N_em: electromagnetic torque (speed dependent)
%   N_friction: friction torque (Coulomb + viscous)

% Friction Model (Wertz Eq. 7-146):
%   N_friction = N_c * sign(s) + f * s
%   N_c: Coulomb friction coefficient
%   f: viscous friction coefficient

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [wheel_state, torque_actual, h_dot] = reaction_wheel_model(wheel_state, torque_cmd, dt, params)

% Varsayılan parametreler
if nargin < 4 || isempty(params)
    params = reaction_wheel_default_params();
end

% Tekerlek durumunu başlat (ilk çağrıda)
if isempty(wheel_state)
    wheel_state = struct();
    wheel_state.speed = 0;  % rpm
    wheel_state.h = 0;      % N·m·s
end

% Mevcut hız
s = wheel_state.speed;  % rpm

% 1. Duty Cycle hesaplama (torque command'dan)
% Komut torku -> duty cycle dönüşümü
% X_dc = torque_cmd / max_torque (saturated)
X_dc = torque_cmd / params.max_torque;
X_dc = max(min(X_dc, 1), -1);  % Saturate to [-1, 1]

% 2. Electromagnetic Torque (Wertz Eq. 7-145)
% N_em(s) ≈ N0 * r * (a - r) / (a * (1 - r))
% where r = 1 - s/s_max for X_dc > 0
%       r = 1 + s/s_max for X_dc < 0
if X_dc >= 0
    r = 1 - s / params.synch_speed;
else
    r = 1 + s / params.synch_speed;
end

% Electromagnetic torque (simplified approximation)
a = params.a_param;  % Wertz parameter
if abs(r) > 0.01 && abs(1 - r) > 0.01
    N_em = params.N0 * r * (a - r) / (a * (1 - r));
else
    N_em = params.N0 * sign(r);  % Near stall or synch speed
end

% Limit N_em to physical bounds
N_em = max(min(N_em, params.max_torque), -params.max_torque);

% Applied electromagnetic torque
N_applied = X_dc * abs(N_em);

% 3. Friction Torque (Wertz Eq. 7-146)
% N_friction = N_c * sign(s) + f * s
if abs(s) > 0.1  % Moving
    N_friction = params.N_c * sign(s) + params.f * s;
else  % Stiction region
    % Static friction (stiction) is higher than kinetic
    N_friction = params.N_c_static * sign(N_applied);
    % Only overcome stiction if applied torque is large enough
    if abs(N_applied) < params.N_c_static
        N_friction = N_applied;  % Stuck, no motion
    end
end

% 4. Net Torque on Wheel
N_net = N_applied - N_friction;

% 5. Quantization (motor driver resolution)
N_net = round(N_net / params.torque_resolution) * params.torque_resolution;

% 6. Torque rate limiting
if isfield(wheel_state, 'last_torque')
    torque_rate = (N_net - wheel_state.last_torque) / dt;
    if abs(torque_rate) > params.max_torque_rate
        N_net = wheel_state.last_torque + sign(torque_rate) * params.max_torque_rate * dt;
    end
end
wheel_state.last_torque = N_net;

% 7. Update wheel speed
% I_w * dω/dt = N_net
% dω = N_net / I_w * dt
d_omega = N_net / params.I_w * dt;  % rad/s
d_speed = d_omega * 60 / (2*pi);    % rpm

wheel_state.speed = s + d_speed;

% Speed saturation
wheel_state.speed = max(min(wheel_state.speed, params.max_speed), -params.max_speed);

% 8. Update angular momentum
% h = I_w * ω
wheel_state.h = params.I_w * wheel_state.speed * 2*pi / 60;  % N·m·s

% Outputs
% Reaction torque on spacecraft = -N_net (Newton's 3rd law)
torque_actual = -N_net;
h_dot = -N_net;  % Rate of change of wheel momentum = -spacecraft torque

end


function params = reaction_wheel_default_params()
% Varsayılan reaction wheel parametreleri
% Tipik küçük uydu reaction wheel değerleri

% Moment of inertia
params.I_w = 0.005;  % kg·m² (typical small RW)

% Synchronous (max) speed
params.synch_speed = 6000;  % rpm
params.max_speed = 6000;    % rpm

% Maximum torque
params.max_torque = 0.02;  % N·m (20 mN·m typical)

% Maximum angular momentum
params.max_momentum = params.I_w * params.max_speed * 2*pi/60;  % N·m·s

% Electromagnetic torque parameters (Wertz Eq. 7-145)
params.N0 = 0.025;   % N·m (max electromagnetic torque)
params.a_param = 0.3; % Wertz 'a' parameter

% Friction coefficients (Wertz Eq. 7-146)
params.N_c = 1e-4;        % N·m Coulomb friction
params.N_c_static = 2e-4; % N·m Static friction (stiction)
params.f = 1e-7;          % N·m/rpm Viscous friction

% Torque resolution (quantization)
params.torque_resolution = 1e-5;  % N·m

% Maximum torque rate
params.max_torque_rate = 0.1;  % N·m/s

% Power consumption (for reference)
params.idle_power = 0.5;   % W
params.max_power = 5.0;    % W

end
