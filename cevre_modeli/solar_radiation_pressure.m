% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Section 5.3.1 - Solar Radiation
% Denklem 5-29

% Bu fonksiyon, güneş radyasyon basıncından kaynaklanan kuvvet ve torku hesaplar.

% r_sat     - uydu pozisyonu ECI frame'de (km) [3x1]
% r_sun     - güneş pozisyonu ECI frame'de (km) [3x1]
% A         - güneşe maruz kalan yüzey alanı (m^2)
% Cr        - reflectivity coefficient (1.0 = tam emici, 2.0 = tam yansıtıcı)
% m         - uydu kütlesi (kg)
% r_cp      - basınç merkezi body frame'de (m) [3x1] (tork hesabı için)
% r_cm      - kütle merkezi body frame'de (m) [3x1] (tork hesabı için)
% in_shadow - gölgede mi? (boolean)
% F_srp     - SRP kuvveti (N) [3x1] - güneş yönünde
% a_srp     - SRP ivmesi (m/s^2) [3x1]
% T_srp     - SRP torku (N·m) [3x1] - basınç merkezi ile kütle merkezi farkından

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: shadow_function

function [F_srp, a_srp, T_srp] = solar_radiation_pressure(r_sat, r_sun, A, Cr, m, r_cp, r_cm)

% Varsayılan değerler
if nargin < 5
    m = 60;             % kg (60kg LEO uydu)
end
if nargin < 6
    r_cp = [0; 0; 0];   % m
end
if nargin < 7
    r_cm = [0; 0; 0];   % m
end

% Sabitler
c = 299792458;          % ışık hızı (m/s)
AU = 149597870.7;       % 1 AU (km)

% Güneş sabitleri (Wertz Denklem 5-29)
F_sun = 1358;           % Ortalama güneş enerji akısı 1 AU'da (W/m^2)

% Gölge kontrolü
[nu, in_shadow] = shadow_function(r_sat, r_sun);

if in_shadow
    % Umbra'da - güneş radyasyonu yok
    F_srp = [0; 0; 0];
    a_srp = [0; 0; 0];
    T_srp = [0; 0; 0];
    return;
end

% Pozisyon vektörlerini sütun vektörüne çevir
if isrow(r_sat)
    r_sat = r_sat';
end
if isrow(r_sun)
    r_sun = r_sun';
end

% Güneş-uydu vektörü
r_sat_sun = r_sun - r_sat;      % km
r_sat_sun_mag = norm(r_sat_sun);
s_hat = r_sat_sun / r_sat_sun_mag;  % Güneş yön birim vektörü

% Güneş akısı (uydu mesafesinde)
% Wertz Denklem 5-29'daki düzeltme terimi (yıllık varyasyon)
% Basitleştirilmiş: 1 AU varsayımı
r_AU = r_sat_sun_mag * 1000 / AU;   % AU cinsinden mesafe (km -> m -> AU)
F_actual = F_sun / r_AU^2;          % W/m^2

% Güneş radyasyon basıncı (Wertz'den)
P_srp = F_actual / c;               % N/m^2 (Pa)

% Penumbra faktörü (0 ile 1 arası)
P_srp = P_srp * nu;

% SRP kuvveti (uydudan güneşe doğru)
% F = -P * Cr * A * s_hat
% Negatif işaret: kuvvet güneşten uzağa doğru
F_srp = -P_srp * Cr * A * s_hat;    % N (ECI frame'de)

% SRP ivmesi
a_srp = F_srp / m;                  % m/s^2

% SRP torku (basınç merkezi ile kütle merkezi arasındaki kol)
r_arm = r_cp - r_cm;                % m (body frame'de)
T_srp = cross(r_arm, F_srp);        % N·m

end
