% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Chapter 4.7 - Effects of Earth's oblateness
% Denklem 4.7.4 (sayfa 178-179)

% Bu fonksiyon, J2 pertürbasyonundan kaynaklanan ivmeyi hesaplar.
% Dünya'nın basıklığı (oblateness) nedeniyle oluşan ek çekimsel ivme.

% R_earth - Dünya'nın ekvator yarıçapı (km)
% J2      - İkinci zonal harmonik katsayısı (boyutsuz)
% mu      - kütleçekimsel parametre (km^3/s^2)
% r       - pozisyon vektörü ECI frame'de (km) [x, y, z]
% a_J2    - J2 pertürbasyon ivmesi ECI frame'de (km/s^2) [ax, ay, az]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function a_J2 = J2_perturbation(r)

global mu

% Sabitler
R_earth = 6378.137;         % Dünya ekvator yarıçapı (km)
J2 = 1.08263e-3;            % İkinci zonal harmonik katsayısı

% Pozisyon bileşenleri
x = r(1);
y = r(2);
z = r(3);

% Pozisyon vektörünün büyüklüğü
r_mag = norm(r);

% J2 pertürbasyon ivmesi bileşenleri (Denklem 4.7.4)
% Ortak çarpan
factor = -3/2 * J2 * mu * R_earth^2 / r_mag^5;

% x bileşeni
ax = factor * x * (1 - 5*(z/r_mag)^2);

% y bileşeni
ay = factor * y * (1 - 5*(z/r_mag)^2);

% z bileşeni
az = factor * z * (3 - 5*(z/r_mag)^2);

% İvme vektörü
a_J2 = [ax; ay; az];

end
