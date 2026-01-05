% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Appendix H - Magnetic Field Models
% Denklem H-1 ile H-14

% Bu fonksiyon, IGRF (International Geomagnetic Reference Field) modelini
% spherical harmonic expansion kullanarak hesaplar.

% r_ECI    - uydu pozisyonu ECI frame'de (km) [3x1] veya [1x3]
% t        - epoch'tan beri geçen süre (saniye)
% nmax     - maksimum derece (varsayılan: 10, IGRF-13 için max 13)
% B_ECI    - manyetik alan vektörü ECI frame'de (nT) [3x1]
% B_NED    - manyetik alan vektörü NED frame'de (nT) [Bn, Be, Bd]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: 
%   get_IGRF_coefficients, legendre_schmidt

function [B_ECI, B_NED] = magnetic_field_IGRF(r_ECI, t, nmax)

% Varsayılan maksimum derece
if nargin < 3
    nmax = 10;
end

% Sabitler
a = 6371.2;                     % IGRF referans yarıçapı (km)
omega_earth = 7.2921159e-5;     % Dünya dönüş hızı (rad/s)

% Pozisyon vektörünü sütun vektörüne çevir
if isrow(r_ECI)
    r_ECI = r_ECI';
end

% ECI'dan spherical koordinatlara dönüşüm
x = r_ECI(1);
y = r_ECI(2);
z = r_ECI(3);

r = norm(r_ECI);
theta = acos(z / r);            % Coelevation (rad)

% Greenwich sidereal time
alpha_G0 = 280.46 * pi/180;     % J2000.0 referans
alpha_G = alpha_G0 + omega_earth * t;

% ECI'dan ECEF'e dönüşüm için longitude
phi_ECI = atan2(y, x);
phi = phi_ECI - alpha_G;        % ECEF longitude (East from Greenwich)
phi = mod(phi, 2*pi);           % 0 ile 2*pi arasında

% IGRF katsayılarını al
[g, h] = get_IGRF_coefficients();

% Manyetik alan bileşenlerini başlat (spherical koordinatlarda)
Br = 0;         % Radyal (dışa pozitif)
Btheta = 0;     % Coelevation (güneye pozitif)
Bphi = 0;       % Azimuthal (doğuya pozitif)

% Spherical harmonic expansion (Wertz Denklem H-12)
for n = 1:nmax
    for m = 0:n
        % Schmidt normalize Legendre fonksiyonları
        [Pnm, dPnm] = legendre_schmidt(n, m, theta);
        
        % Gauss katsayıları
        gnm = g(n, m+1);
        hnm = h(n, m+1);
        
        % Trigonometrik terimler
        cos_mphi = cos(m * phi);
        sin_mphi = sin(m * phi);
        
        % Ortak çarpan
        factor = (a/r)^(n+2);
        
        % Alan bileşenleri (Wertz Denklem H-12)
        Br = Br + (n+1) * factor * (gnm*cos_mphi + hnm*sin_mphi) * Pnm;
        Btheta = Btheta - factor * (gnm*cos_mphi + hnm*sin_mphi) * dPnm;
        Bphi = Bphi + factor * m * (-gnm*sin_mphi + hnm*cos_mphi) * Pnm / sin(theta);
    end
end

% NED koordinatlarına dönüşüm (Wertz Denklem H-13)
% Basitleştirilmiş: geodetic = geocentric varsayımı
Bn = -Btheta;           % Kuzey
Be = Bphi;              % Doğu
Bd = -Br;               % Aşağı (nadir)

B_NED = [Bn; Be; Bd];

% ECI koordinatlarına dönüşüm (Wertz Denklem H-14)
cos_theta = cos(theta);
sin_theta = sin(theta);
cos_phi_eci = cos(phi_ECI);
sin_phi_eci = sin(phi_ECI);

Bx = (Br*cos_theta + Btheta*sin_theta)*cos_phi_eci - Bphi*sin_phi_eci;
By = (Br*cos_theta + Btheta*sin_theta)*sin_phi_eci + Bphi*cos_phi_eci;
Bz = Br*sin_theta - Btheta*cos_theta;

B_ECI = [Bx; By; Bz];

end


function [g, h] = get_IGRF_coefficients()
% IGRF-13 katsayıları (2020 epoch, nT)
% Kaynak: https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field
% Sadece n=1:10 için (yeterli doğruluk)

% g(n, m+1) ve h(n, m+1) formatında
% Not: m=0 için h=0

g = zeros(13, 14);
h = zeros(13, 14);

% n=1
g(1,1) = -29404.8; g(1,2) = -1450.9;
h(1,1) = 0;        h(1,2) = 4652.5;

% n=2
g(2,1) = -2499.6; g(2,2) = 2982.0; g(2,3) = 1677.0;
h(2,1) = 0;       h(2,2) = -2991.6; h(2,3) = -734.6;

% n=3
g(3,1) = 1363.2; g(3,2) = -2381.2; g(3,3) = 1236.2; g(3,4) = 525.7;
h(3,1) = 0;      h(3,2) = -82.1;   h(3,3) = 241.9;  h(3,4) = -543.4;

% n=4
g(4,1) = 903.0; g(4,2) = 809.5; g(4,3) = 86.3; g(4,4) = -309.4; g(4,5) = 48.0;
h(4,1) = 0;     h(4,2) = 281.9; h(4,3) = -158.4; h(4,4) = 199.7; h(4,5) = -349.7;

% n=5
g(5,1) = -234.3; g(5,2) = 363.2; g(5,3) = 47.7; g(5,4) = 187.8; g(5,5) = -140.7; g(5,6) = -151.2;
h(5,1) = 0;      h(5,2) = 46.9;  h(5,3) = 196.5; h(5,4) = -119.3; h(5,5) = 16.0;  h(5,6) = 100.1;

% n=6
g(6,1) = 66.0; g(6,2) = 65.5; g(6,3) = -19.1; g(6,4) = 72.9; g(6,5) = -121.5; g(6,6) = -36.2; g(6,7) = 13.5;
h(6,1) = 0;    h(6,2) = -20.8; h(6,3) = 59.9; h(6,4) = -82.6; h(6,5) = -27.2; h(6,6) = 6.4;   h(6,7) = -53.1;

% n=7
g(7,1) = 72.9; g(7,2) = -54.1; g(7,3) = -6.0; g(7,4) = -19.5; g(7,5) = 16.5; g(7,6) = 3.9; g(7,7) = -24.3; g(7,8) = -3.2;
h(7,1) = 0;    h(7,2) = 24.2;  h(7,3) = -4.0; h(7,4) = 2.8;   h(7,5) = 21.0; h(7,6) = -6.2; h(7,7) = 10.2;  h(7,8) = -20.3;

% n=8
g(8,1) = 3.2; g(8,2) = -5.9; g(8,3) = 3.6; g(8,4) = -5.0; g(8,5) = -0.8; g(8,6) = 0.9; g(8,7) = -8.0; g(8,8) = -1.1; g(8,9) = 4.7;
h(8,1) = 0;   h(8,2) = 8.4;  h(8,3) = -1.0; h(8,4) = -14.0; h(8,5) = 9.4; h(8,6) = 5.0; h(8,7) = 3.3; h(8,8) = -9.5; h(8,9) = -0.5;

% n=9
g(9,1) = 5.3; g(9,2) = 3.9; g(9,3) = -8.0; g(9,4) = -0.1; g(9,5) = -2.3; g(9,6) = -5.6; g(9,7) = 2.9; g(9,8) = 3.3; g(9,9) = 0.1; g(9,10) = -1.5;
h(9,1) = 0;   h(9,2) = -2.0; h(9,3) = 0.5; h(9,4) = 4.9; h(9,5) = -3.0; h(9,6) = -0.5; h(9,7) = 0.0; h(9,8) = -2.6; h(9,9) = 1.9; h(9,10) = -0.4;

% n=10
g(10,1) = -2.1; g(10,2) = -0.5; g(10,3) = 0.6; g(10,4) = 1.3; g(10,5) = -1.6; g(10,6) = -0.9; g(10,7) = 0.7; g(10,8) = 0.3; g(10,9) = 1.7; g(10,10) = -0.2; g(10,11) = 0.7;
h(10,1) = 0;   h(10,2) = 2.9; h(10,3) = -2.1; h(10,4) = -0.4; h(10,5) = 0.2; h(10,6) = 0.2; h(10,7) = 0.6; h(10,8) = -0.8; h(10,9) = 0.5; h(10,10) = 1.4; h(10,11) = -1.2;

end


function [Pnm, dPnm] = legendre_schmidt(n, m, theta)
% Schmidt quasi-normalized associated Legendre functions
% Wertz Denklem H-4 ve H-10

cos_theta = cos(theta);
sin_theta = sin(theta);

% Recursion için başlangıç değerleri
if n == 0 && m == 0
    Pnm = 1;
    dPnm = 0;
    return;
end

% P ve dP hesaplama (recursion ile)
P = zeros(n+1, n+1);
P(1,1) = 1;

for i = 1:n
    % Diagonal terms: P(i,i)
    P(i+1, i+1) = sin_theta * P(i, i) * sqrt((2*i-1)/(2*i));
    
    % Sub-diagonal terms: P(i,i-1)
    if i > 1
        P(i+1, i) = cos_theta * P(i, i) * sqrt(2*i-1);
    end
end

% Diğer terimler
for i = 2:n
    for j = 0:i-2
        K = ((i-1)^2 - j^2) / ((2*i-1)*(2*i-3));
        if i == 1
            P(i+1, j+1) = cos_theta * P(i, j+1);
        else
            P(i+1, j+1) = cos_theta * P(i, j+1) - K * P(i-1, j+1);
        end
    end
end

Pnm = P(n+1, m+1);

% Türev hesaplama (Wertz Denklem H-10)
if m == n
    if n == 0
        dPnm = 0;
    else
        dPnm = sin_theta * legendre_schmidt_value(n-1, n-1, theta) + cos_theta * P(n, n);
    end
else
    K = ((n-1)^2 - m^2) / ((2*n-1)*(2*n-3));
    if n >= 2
        dPnm = cos_theta * legendre_schmidt_derivative(n-1, m, theta) ...
               - sin_theta * P(n, m+1) - K * legendre_schmidt_derivative(n-2, m, theta);
    else
        dPnm = cos_theta * legendre_schmidt_derivative(n-1, m, theta) - sin_theta * P(n, m+1);
    end
end

end


function Pnm = legendre_schmidt_value(n, m, theta)
% Yardımcı fonksiyon: sadece değer döndür
[Pnm, ~] = legendre_schmidt(n, m, theta);
end


function dPnm = legendre_schmidt_derivative(n, m, theta)
% Yardımcı fonksiyon: sadece türev döndür
if n < 0
    dPnm = 0;
    return;
end
[~, dPnm] = legendre_schmidt(n, m, theta);
end
