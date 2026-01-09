% Sun Position Calculation from UTC Time
% Astronomik almanac formülleri kullanılarak güneş konumu hesaplanır.
%
% Bu fonksiyon, GNSS'den gelen UTC zaman bilgisiyle güneşin ECI 
% koordinatlarındaki konumunu hesaplar.
%
% Referanslar:
%   - Vallado (2013), "Fundamentals of Astrodynamics and Applications"
%   - Meeus (1998), "Astronomical Algorithms"
%   - USNO/SOFA Almanac

% Girdiler:
%   utc_time  - UTC zaman [year, month, day, hour, min, sec] veya datetime
%   r_sat     - Uydu pozisyonu ECI (m) [3x1] (opsiyonel, eclipse için)
%
% Çıktılar:
%   s_eci     - Güneş birim vektörü ECI frame'de [3x1]
%   s_eci_km  - Güneş pozisyonu ECI (km) [3x1]
%   eclipse   - Eclipse durumu (0: güneşte, 1: gölgede)
%   sun_dist  - Güneş-Dünya mesafesi (AU)

% Güneş konumu hesaplama adımları:
%   1. Julian Date hesapla
%   2. Mean anomaly, ecliptic longitude hesapla
%   3. Ecliptic -> ECI dönüşümü
%   4. Eclipse kontrolü (opsiyonel)

function [s_eci, s_eci_km, eclipse, sun_dist] = sun_position_eci(utc_time, r_sat)

% Giriş kontrolü
if isa(utc_time, 'datetime')
    year = utc_time.Year;
    month = utc_time.Month;
    day = utc_time.Day;
    hour = utc_time.Hour;
    minute = utc_time.Minute;
    second = utc_time.Second;
elseif length(utc_time) >= 6
    year = utc_time(1);
    month = utc_time(2);
    day = utc_time(3);
    hour = utc_time(4);
    minute = utc_time(5);
    second = utc_time(6);
else
    error('utc_time must be [year, month, day, hour, min, sec] or datetime');
end

% Eclipse kontrolü için uydu pozisyonu
if nargin < 2
    r_sat = [];
end

%% ==================== JULIAN DATE ====================
% Julian Date hesaplama (Vallado Eq. 3-13)
JD = julian_date(year, month, day, hour, minute, second);

% Julian centuries from J2000.0
T_UT1 = (JD - 2451545.0) / 36525;

%% ==================== SUN POSITION (Low Precision) ====================
% Vallado Algorithm 29 (pg 279-280) - Accuracy ~0.01 deg

% Mean longitude of Sun (deg)
lambda_M_sun = 280.4606184 + 36000.77005361 * T_UT1;
lambda_M_sun = mod(lambda_M_sun, 360);

% Mean anomaly of Sun (deg)
M_sun = 357.5277233 + 35999.05034 * T_UT1;
M_sun = mod(M_sun, 360);
M_sun_rad = deg2rad(M_sun);

% Ecliptic longitude of Sun (deg)
lambda_ecliptic = lambda_M_sun + 1.914666471 * sin(M_sun_rad) + ...
                  0.019994643 * sin(2 * M_sun_rad);
lambda_ecliptic = mod(lambda_ecliptic, 360);
lambda_ecliptic_rad = deg2rad(lambda_ecliptic);

% Distance to Sun (AU)
sun_dist = 1.000140612 - 0.016708617 * cos(M_sun_rad) - ...
           0.000139589 * cos(2 * M_sun_rad);

% Obliquity of ecliptic (deg)
epsilon = 23.439291 - 0.0130042 * T_UT1;
epsilon_rad = deg2rad(epsilon);

%% ==================== ECLIPTIC TO ECI ====================
% Sun unit vector in ECI (Vallado Eq. 3-87)
s_eci = [cos(lambda_ecliptic_rad);
         cos(epsilon_rad) * sin(lambda_ecliptic_rad);
         sin(epsilon_rad) * sin(lambda_ecliptic_rad)];

% Normalize (should already be ~1)
s_eci = s_eci / norm(s_eci);

% Sun position in km (1 AU = 149597870.7 km)
AU_to_km = 149597870.7;
s_eci_km = s_eci * sun_dist * AU_to_km;

%% ==================== ECLIPSE CHECK ====================
eclipse = 0;  % Default: not in eclipse

if ~isempty(r_sat)
    % Uydu pozisyonunu km'ye çevir (m olarak geldiyse)
    if norm(r_sat) > 1e6  % Muhtemelen metre
        r_sat_km = r_sat / 1000;
    else
        r_sat_km = r_sat;
    end
    
    eclipse = check_eclipse(r_sat_km, s_eci_km);
end

end


%% ==================== HELPER FUNCTIONS ====================

function JD = julian_date(year, month, day, hour, minute, second)
% Julian Date hesaplama (Vallado Eq. 3-13)
    
    % Fractional day
    day_frac = day + (hour + minute/60 + second/3600) / 24;
    
    % Julian Date formula
    JD = 367 * year - floor(7 * (year + floor((month + 9) / 12)) / 4) + ...
         floor(275 * month / 9) + day_frac + 1721013.5;
end


function eclipse = check_eclipse(r_sat_km, s_sun_km)
% Cylindrical shadow model (Vallado pg 301)
% Dünya'nın gölgesinde olup olmadığını kontrol eder
    
    % Earth radius (km)
    R_earth = 6378.137;
    
    % Sun direction unit vector
    s_hat = s_sun_km / norm(s_sun_km);
    
    % Satellite position projection onto sun direction
    sat_sun_proj = dot(r_sat_km, s_hat);
    
    % Satellite is on the sun side of Earth
    if sat_sun_proj > 0
        eclipse = 0;
        return;
    end
    
    % Perpendicular distance from satellite to sun-earth line
    r_perp = norm(r_sat_km - sat_sun_proj * s_hat);
    
    % Simple cylindrical shadow model
    if r_perp < R_earth && sat_sun_proj < 0
        eclipse = 1;  % In Earth's shadow
    else
        eclipse = 0;  % In sunlight
    end
end


%% ==================== ADDITIONAL UTILITY FUNCTIONS ====================

function s_body = sun_vector_body(s_eci, q)
% Güneş vektörünü body frame'e dönüştür
% Bu fonksiyon, kestirilen quaternion ile güneşin body frame'deki
% konumunu hesaplar.
%
% Girdiler:
%   s_eci - Güneş birim vektörü ECI [3x1]
%   q     - Attitude quaternion (ECI -> Body) [4x1]
%
% Çıktı:
%   s_body - Güneş birim vektörü body frame'de [3x1]

    % Quaternion'dan DCM
    A = quaternion_to_dcm_sun(q);
    
    % ECI -> Body dönüşümü
    s_body = A * s_eci;
end


function A = quaternion_to_dcm_sun(q)
% Quaternion'dan DCM'e dönüşüm
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    A = [q1^2-q2^2-q3^2+q4^2,   2*(q1*q2+q3*q4),       2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4),       -q1^2+q2^2-q3^2+q4^2,  2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4),       2*(q2*q3-q1*q4),       -q1^2-q2^2+q3^2+q4^2];
end
