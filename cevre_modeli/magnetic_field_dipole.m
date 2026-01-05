% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Appendix H - Magnetic Field Models
% Denklem H-17 ile H-25

% Bu fonksiyon, Dünya'nın manyetik alanını dipol modeli kullanarak hesaplar.
% ECI (Geocentric Inertial) koordinatlarında B vektörünü döndürür.

% r_ECI    - uydu pozisyonu ECI frame'de (km) [3x1] veya [1x3]
% t        - epoch'tan beri geçen süre (saniye)
% B_ECI    - manyetik alan vektörü ECI frame'de (nT) [3x1]
% B_mag    - manyetik alan şiddeti (nT)

% Dipol parametreleri (IGRF-13, 2020 epoch için güncellenmiş):
% g10 = -29404.8 nT
% g11 = -1450.9 nT  
% h11 = 4652.5 nT
% Dipol momenti: H0 = sqrt(g10^2 + g11^2 + h11^2) ≈ 29873 nT (yüzeyde)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [B_ECI, B_mag] = magnetic_field_dipole(r_ECI, t)

% Sabitler
a = 6371.2;                     % Dünya ortalama yarıçapı (km) - IGRF referans

% IGRF-13 (2020 epoch) Gauss katsayıları (nT)
% Kaynak: https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field
g10 = -29404.8;
g11 = -1450.9;
h11 = 4652.5;

% Dipol momenti büyüklüğü (Wertz Denklem H-18)
H0 = sqrt(g10^2 + g11^2 + h11^2);   % nT

% Dipol yönü (Wertz Denklem H-19 ve H-20)
theta_m = acos(g10 / H0);           % Dipol coelevation (rad)
phi_m = atan2(h11, g11);            % Dipol East longitude (rad)

% Greenwich sidereal time hesaplama
% Referans: J2000.0 epoch (1 Ocak 2000, 12:00 UT)
omega_earth = 7.2921159e-5;         % Dünya dönüş hızı (rad/s)

% alpha_G0: Greenwich'in sağ açıklığı referans zamanında
% J2000.0 için: alpha_G0 ≈ 280.46 derece
alpha_G0 = 280.46 * pi/180;         % rad

% Greenwich sidereal time
alpha_G = alpha_G0 + omega_earth * t;   % rad

% Dipol birim vektörü ECI'da (Wertz Denklem H-23)
alpha_m = alpha_G + phi_m;          % Dipol sağ açıklığı
m_hat = [sin(theta_m) * cos(alpha_m);
         sin(theta_m) * sin(alpha_m);
         cos(theta_m)];

% Pozisyon vektörünü sütun vektörüne çevir
if isrow(r_ECI)
    r_ECI = r_ECI';
end

% Pozisyon büyüklüğü ve birim vektör
R = norm(r_ECI);
R_hat = r_ECI / R;

% Manyetik alan hesaplama (Wertz Denklem H-22)
% B = (a^3 * H0 / R^3) * [3*(m_hat · R_hat)*R_hat - m_hat]
m_dot_R = dot(m_hat, R_hat);
B_ECI = (a^3 * H0 / R^3) * (3 * m_dot_R * R_hat - m_hat);

% Manyetik alan şiddeti
B_mag = norm(B_ECI);

end
