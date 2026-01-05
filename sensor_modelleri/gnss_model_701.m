% GNSS-701 Satellite GNSS Receiver Model
% AAC Clyde Space GNSS-701 Datasheet parametreleri kullanılmıştır.

% Bu fonksiyon, GNSS-701 GPS alıcı ölçümlerini simüle eder.

% r_true     - gerçek pozisyon ECI frame'de (km) [3x1]
% v_true     - gerçek hız ECI frame'de (km/s) [3x1]
% t_true     - gerçek zaman (s)
% params     - sensör parametreleri (opsiyonel, varsayılan GNSS-701)
% gnss_state - sensör durumu (clock bias için, opsiyonel)
% r_meas     - ölçülen pozisyon (km) [3x1]
% v_meas     - ölçülen hız (km/s) [3x1]
% t_meas     - ölçülen zaman (s)
% valid      - ölçüm geçerliliği (boolean)

% GNSS-701 Single Frequency Specifications (Datasheet):
%   Position Accuracy: 1.5 meters
%   Velocity Accuracy: <0.03 m/s RMS
%   Time Accuracy: 20 ns RMS
%   Frequencies: L1, E1, B1
%   Constellations: GPS, GPS+GLO, GPS+GLO+GAL
%   Time to First Fix: Cold Start 40s, Hot Start 19s
%   Update Rate: 1 Hz (PPS output)

% GNSS ölçüm modeli:
%   r_meas = r_true + n_pos
%   v_meas = v_true + n_vel
%   t_meas = t_true + bias_clock + n_time

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [r_meas, v_meas, t_meas, valid, gnss_state] = gnss_model_701(r_true, v_true, t_true, params, gnss_state)

% Varsayılan parametreler (GNSS-701 datasheet'ten)
if nargin < 4 || isempty(params)
    params = gnss701_default_params();
end

% GNSS durumunu başlat (ilk çağrıda)
if nargin < 5 || isempty(gnss_state)
    gnss_state = struct();
    gnss_state.clock_bias = params.time_accuracy * randn;  % Başlangıç saat bias
    gnss_state.clock_drift = 1e-9 * randn;  % s/s clock drift
    gnss_state.last_update = t_true;
    gnss_state.fix_acquired = false;
    gnss_state.fix_time = 0;
end

% Vektörleri sütun vektörüne çevir
if isrow(r_true)
    r_true = r_true';
end
if isrow(v_true)
    v_true = v_true';
end

% Update rate kontrolü (1 Hz)
dt = t_true - gnss_state.last_update;
if dt < (1/params.update_rate) && gnss_state.fix_acquired
    % Henüz güncelleme zamanı değil, son değerleri kullan
    r_meas = gnss_state.last_r;
    v_meas = gnss_state.last_v;
    t_meas = gnss_state.last_t;
    valid = gnss_state.fix_acquired;
    return;
end

% Fix kontrolü (TTFF - Time to First Fix)
if ~gnss_state.fix_acquired
    gnss_state.fix_time = gnss_state.fix_time + dt;
    if gnss_state.fix_time < params.ttff_cold
        % Henüz fix alınmadı
        r_meas = NaN(3,1);
        v_meas = NaN(3,1);
        t_meas = NaN;
        valid = false;
        gnss_state.last_update = t_true;
        return;
    else
        gnss_state.fix_acquired = true;
    end
end

% Clock bias update (random walk)
gnss_state.clock_bias = gnss_state.clock_bias + ...
                        gnss_state.clock_drift * dt + ...
                        params.clock_noise * sqrt(dt) * randn;

% 1. Position Measurement
% Position Accuracy: 1.5 m = 0.0015 km (1-sigma, 3D)
% 3D accuracy → each axis: 1.5/sqrt(3) m
pos_sigma = params.position_accuracy / sqrt(3);  % km, per axis
pos_noise = pos_sigma * randn(3,1);
r_meas = r_true + pos_noise;

% 2. Velocity Measurement
% Velocity Accuracy: 0.03 m/s = 0.00003 km/s (1-sigma, 3D)
vel_sigma = params.velocity_accuracy / sqrt(3);  % km/s, per axis
vel_noise = vel_sigma * randn(3,1);
v_meas = v_true + vel_noise;

% 3. Time Measurement
% Time Accuracy: 20 ns RMS
time_noise = params.time_accuracy * randn;
t_meas = t_true + gnss_state.clock_bias + time_noise;

% Geçerlilik kontrolü
% Yükseklik kontrolü (LEO için tipik olarak 200-2000 km)
r_mag = norm(r_true);
altitude = r_mag - 6378.137;  % km (WGS84 equatorial radius)
valid = (altitude > 200) && (altitude < 2000);

% LEO'da GPS görünürlüğü bazen sınırlı olabilir
% Basit model: %95 availability
if rand > params.availability
    valid = false;
    r_meas = NaN(3,1);
    v_meas = NaN(3,1);
    t_meas = NaN;
end

% State güncelle
gnss_state.last_update = t_true;
gnss_state.last_r = r_meas;
gnss_state.last_v = v_meas;
gnss_state.last_t = t_meas;

end


function params = gnss701_default_params()
% GNSS-701 varsayılan parametreleri

% Position Accuracy: 1.5 m
params.position_accuracy = 1.5 / 1000;  % km

% Velocity Accuracy: 0.03 m/s
params.velocity_accuracy = 0.03 / 1000;  % km/s

% Time Accuracy: 20 ns
params.time_accuracy = 20e-9;  % s

% Time to First Fix: Cold Start 40s, Hot Start 19s
params.ttff_cold = 40;  % s
params.ttff_hot = 19;   % s

% Update Rate: 1 Hz
params.update_rate = 1;  % Hz

% Clock noise (Allan variance based estimate)
params.clock_noise = 1e-9;  % s/√s

% Availability (LEO'da tipik GPS availability)
params.availability = 0.95;  % 95%

% Dilution of Precision (DOP) - tipik değer
params.gdop_typical = 2.0;

% Maximum altitude for operation
params.max_altitude = 3000;  % km

end
