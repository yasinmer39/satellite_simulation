% Orbital Mechanics for Engineering Students (Curtis) kitabından alınmıştır.
% Section 10.5-10.6 - Atmospheric Drag

% Bu fonksiyon, atmosferik sürüklenme (drag) ivmesini hesaplar.

% r_sat     - uydu pozisyonu ECI frame'de (km) [3x1]
% v_sat     - uydu hızı ECI frame'de (km/s) [3x1]
% Cd        - drag katsayısı (tipik: 2.0-2.5)
% A         - çapraz kesit alanı (m^2)
% m         - uydu kütlesi (kg)
% a_drag    - drag ivmesi ECI frame'de (km/s^2) [3x1]
% rho       - atmosfer yoğunluğu (kg/m^3)

% Atmosfer modeli: Üstel (exponential) model
% Curtis Tablo 10.1'den interpolasyon

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: atmospheric_density

function [a_drag, rho] = atmospheric_drag(r_sat, v_sat, Cd, A, m)

% Varsayılan değerler
if nargin < 3
    Cd = 2.2;           % Tipik değer
end
if nargin < 4
    A = 1.0;            % m^2 (varsayılan çapraz kesit)
end
if nargin < 5
    m = 60;             % kg (60kg LEO uydu)
end

% Sabitler
R_earth = 6378.137;     % Dünya yarıçapı (km)
omega_earth = 7.2921159e-5;  % Dünya dönüş hızı (rad/s)

% Pozisyon vektörlerini sütun vektörüne çevir
if isrow(r_sat)
    r_sat = r_sat';
end
if isrow(v_sat)
    v_sat = v_sat';
end

% Yükseklik hesapla
r_mag = norm(r_sat);
h = r_mag - R_earth;    % km

% Yükseklik kontrolü (atmosfer dışı)
if h > 1000
    a_drag = [0; 0; 0];
    rho = 0;
    return;
end

% Atmosfer yoğunluğu
rho = atmospheric_density(h);   % kg/m^3

% Atmosfere göre bağıl hız
% Atmosfer Dünya ile birlikte döner
% v_atm = omega_earth × r_sat
omega_vec = [0; 0; omega_earth];    % rad/s
v_atm = cross(omega_vec, r_sat);    % km/s

% Bağıl hız
v_rel = v_sat - v_atm;              % km/s
v_rel_mag = norm(v_rel);            % km/s
v_rel_m = v_rel_mag * 1000;         % m/s

% Drag ivmesi (Curtis Denklem 10.32)
% a_drag = -½ * rho * v² * Cd * A / m * v_hat
% Birim: km/s²
a_drag_mag = 0.5 * rho * v_rel_m^2 * Cd * A / m;   % m/s^2
a_drag_mag = a_drag_mag / 1000;     % km/s^2

% Yön: hıza karşı
v_rel_hat = v_rel / v_rel_mag;
a_drag = -a_drag_mag * v_rel_hat;   % km/s^2

end


function rho = atmospheric_density(h)
% Üstel atmosfer modeli
% Curtis Tablo 10.1'den

% Yükseklik bantları ve parametreler
% [h0 (km), rho0 (kg/m³), H (km)]
atm_data = [
    0,      1.225,      7.249;
    25,     3.899e-2,   6.349;
    30,     1.774e-2,   6.682;
    40,     3.972e-3,   7.554;
    50,     1.057e-3,   8.382;
    60,     3.206e-4,   7.714;
    70,     8.770e-5,   6.549;
    80,     1.905e-5,   5.799;
    90,     3.396e-6,   5.382;
    100,    5.297e-7,   5.877;
    110,    9.661e-8,   7.263;
    120,    2.438e-8,   9.473;
    130,    8.484e-9,   12.636;
    140,    3.845e-9,   16.149;
    150,    2.070e-9,   22.523;
    180,    5.464e-10,  29.740;
    200,    2.789e-10,  37.105;
    250,    7.248e-11,  45.546;
    300,    2.418e-11,  53.628;
    350,    9.518e-12,  53.298;
    400,    3.725e-12,  58.515;
    450,    1.585e-12,  60.828;
    500,    6.967e-13,  63.822;
    600,    1.454e-13,  71.835;
    700,    3.614e-14,  88.667;
    800,    1.170e-14,  124.64;
    900,    5.245e-15,  181.05;
    1000,   3.019e-15,  268.00;
];

% Yükseklik sınır kontrolü
if h < 0
    h = 0;
end
if h > 1000
    rho = 0;
    return;
end

% Uygun bantı bul
n = size(atm_data, 1);
for i = 1:n-1
    if h >= atm_data(i, 1) && h < atm_data(i+1, 1)
        h0 = atm_data(i, 1);
        rho0 = atm_data(i, 2);
        H = atm_data(i, 3);
        break;
    end
end

% Son bant için
if h >= atm_data(n, 1)
    h0 = atm_data(n, 1);
    rho0 = atm_data(n, 2);
    H = atm_data(n, 3);
end

% Üstel formül
rho = rho0 * exp(-(h - h0) / H);

end
