% Gas Jet (Thruster) Actuator Model
% Wertz Chapter 7.10 - Modeling Gas-Jet Control Systems referans alınmıştır.

% Bu fonksiyon, gas jet (cold gas / hot gas thruster) davranışını simüle eder.
% Thrust profile, delays, rise/fall times modellenir.

% inputs:
%   cmd_on      - açma komutu (boolean)
%   t           - mevcut zaman (s)
%   jet_state   - thruster durumu struct
%   params      - thruster parametreleri (opsiyonel)
%
% outputs:
%   jet_state   - güncellenmiş thruster durumu
%   thrust      - anlık itme kuvveti (N)
%   torque      - uzay aracına uygulanan tork (N·m) [3x1]

% Thrust Profile (Wertz Fig. 7-30):
%   t0: commanded start time
%   t1: thrust begins buildup
%   t2: thrust reaches steady state
%   t3: commanded stop time
%   t4: thrust begins decay
%   t5: thrust reaches zero
%
%   Delays: t1-t0, t4-t3 (valve delays)
%   Rise time: t2-t1
%   Fall time: t5-t4

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [jet_state, thrust, torque] = gas_jet_model(cmd_on, t, jet_state, params)

% Varsayılan parametreler
if nargin < 4 || isempty(params)
    params = gas_jet_default_params();
end

% Jet durumunu başlat (ilk çağrıda)
if isempty(jet_state)
    jet_state = struct();
    jet_state.cmd_on = false;
    jet_state.firing = false;
    jet_state.thrust_level = 0;
    jet_state.cmd_on_time = -inf;
    jet_state.cmd_off_time = -inf;
    jet_state.total_impulse = 0;
    jet_state.fuel_consumed = 0;
end

% Komut değişikliği algılama
if cmd_on && ~jet_state.cmd_on
    % Rising edge - açma komutu
    jet_state.cmd_on_time = t;
    jet_state.cmd_on = true;
elseif ~cmd_on && jet_state.cmd_on
    % Falling edge - kapama komutu
    jet_state.cmd_off_time = t;
    jet_state.cmd_on = false;
end

% Thrust profile hesaplama
if jet_state.cmd_on
    % Açık durumda
    t_since_cmd = t - jet_state.cmd_on_time;
    
    if t_since_cmd < params.delay_on
        % Henüz valve açılmadı
        thrust_fraction = 0;
    elseif t_since_cmd < params.delay_on + params.rise_time
        % Rise phase (linear ramp)
        t_rise = t_since_cmd - params.delay_on;
        thrust_fraction = t_rise / params.rise_time;
    else
        % Steady state
        thrust_fraction = 1.0;
    end
else
    % Kapalı durumda
    t_since_cmd = t - jet_state.cmd_off_time;
    
    if t_since_cmd < params.delay_off
        % Henüz valve kapanmadı, hala steady state
        thrust_fraction = 1.0;
    elseif t_since_cmd < params.delay_off + params.fall_time
        % Fall phase (linear decay)
        t_fall = t_since_cmd - params.delay_off;
        thrust_fraction = 1 - t_fall / params.fall_time;
    else
        % Tamamen kapalı
        thrust_fraction = 0;
    end
end

% Thrust hesaplama
thrust_fraction = max(min(thrust_fraction, 1), 0);

% Minimum impulse bit kontrolü
if thrust_fraction > 0 && thrust_fraction < params.min_thrust_fraction
    thrust_fraction = params.min_thrust_fraction;
end

% Nominal thrust with noise
thrust_noise = 1 + params.thrust_noise * randn;
thrust = params.nominal_thrust * thrust_fraction * thrust_noise;

% Thrust quantization (PWM etkisi)
thrust = round(thrust / params.thrust_resolution) * params.thrust_resolution;

% Torque hesaplama
% N = r × F
% r: thruster pozisyonu body frame'de
% F: thrust vektörü body frame'de
force_vector = params.thrust_direction * thrust;
torque = cross(params.position, force_vector);

% State güncelle
jet_state.thrust_level = thrust;
jet_state.firing = (thrust > 0);

% Yakıt tüketimi ve total impulse
if jet_state.firing
    dt = 0.001;  % Varsayılan dt (gerçekte dışarıdan gelmeli)
    jet_state.total_impulse = jet_state.total_impulse + thrust * dt;
    jet_state.fuel_consumed = jet_state.fuel_consumed + thrust * dt / params.Isp / 9.80665;
end

end


function params = gas_jet_default_params()
% Varsayılan gas jet parametreleri
% Cold gas thruster (tipik küçük uydu)

% Nominal thrust
params.nominal_thrust = 1.0;  % N

% Specific impulse
params.Isp = 70;  % s (cold gas, N2 or similar)

% Timing parameters (Wertz Fig. 7-30)
params.delay_on = 0.010;   % s (valve opening delay)
params.delay_off = 0.005;  % s (valve closing delay)
params.rise_time = 0.020;  % s
params.fall_time = 0.015;  % s

% Minimum impulse bit
params.min_impulse_bit = 0.001;  % N·s
params.min_thrust_fraction = 0.1;

% Thrust noise (1-sigma, fraction)
params.thrust_noise = 0.02;  % 2%

% Thrust resolution
params.thrust_resolution = 0.01;  % N

% Thruster geometry (body frame)
params.position = [0.3; 0; 0];         % m (thruster position from CoM)
params.thrust_direction = [0; 1; 0];   % unit vector (thrust direction)

% Fuel parameters
params.initial_fuel_mass = 1.0;  % kg
params.tank_pressure_initial = 2e6;  % Pa (20 bar)

end
