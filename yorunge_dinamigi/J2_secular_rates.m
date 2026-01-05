% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Chapter 4.7 - Effects of Earth's oblateness
% Denklem 4.7.13 ve 4.7.14 (sayfa 185-186)

% Bu fonksiyon, J2 pertürbasyonunun RAAN ve argument of perigee üzerindeki
% seküler (uzun vadeli) etkilerini hesaplar.

% a         - yarı büyük eksen (km)
% e         - eksantriklik
% incl      - eğilim (rad)
% RAAN_dot  - RAAN değişim hızı (rad/s)
% w_dot     - argument of perigee değişim hızı (rad/s)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [RAAN_dot, w_dot] = J2_secular_rates(a, e, incl)

global mu

% Sabitler
R_earth = 6378.137;         % Dünya ekvator yarıçapı (km)
J2 = 1.08263e-3;            % İkinci zonal harmonik katsayısı

% Mean motion (ortalama açısal hız)
n = sqrt(mu / a^3);         % rad/s

% Semi-latus rectum
p = a * (1 - e^2);          % km

% RAAN değişim hızı (Denklem 4.7.13)
% Negatif işaret: RAAN batıya doğru kayar (prograd yörüngeler için)
RAAN_dot = -3/2 * J2 * (R_earth/p)^2 * n * cos(incl);   % rad/s

% Argument of perigee değişim hızı (Denklem 4.7.14)
w_dot = 3/4 * J2 * (R_earth/p)^2 * n * (5*cos(incl)^2 - 1);  % rad/s

end
