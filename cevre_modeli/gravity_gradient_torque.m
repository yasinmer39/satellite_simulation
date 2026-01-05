% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Section 5.2 - The Earth's Gravitational Field
% Denklem 5-26

% Bu fonksiyon, kütleçekimsel gradyent torkunu hesaplar.
% Bu tork, uydunun atalet momentlerinin eşit olmamasından kaynaklanır.

% r_sat     - uydu pozisyonu ECI frame'de (km) [3x1]
% I         - atalet matrisi body frame'de (kg·m²) [3x3]
% A_bn      - body frame'den ECI frame'e dönüşüm matrisi [3x3]
%             (alternatif: quaternion veya Euler açıları)
% T_gg      - gravity gradient torku body frame'de (N·m) [3x1]

% Not: Bu tork uydunun en büyük atalet eksenini nadir'e (Dünya merkezine)
%      doğru çevirmeye çalışır.

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function T_gg = gravity_gradient_torque(r_sat, I, A_bn)

global mu

% mu tanımlı değilse varsayılan değer
if isempty(mu)
    mu = 398600;        % km³/s²
end

% Pozisyon vektörünü sütun vektörüne çevir
if isrow(r_sat)
    r_sat = r_sat';
end

% Pozisyon büyüklüğü
r_mag = norm(r_sat);

% Nadir birim vektörü ECI frame'de (Dünya'ya doğru)
n_ECI = -r_sat / r_mag;

% Nadir vektörünü body frame'e dönüştür
% n_body = A_bn' * n_ECI (ECI'dan body'ye)
n_body = A_bn' * n_ECI;

% Gravity gradient torku (Wertz Denklem 5-26)
% T_gg = (3*mu/r³) * n × (I·n)
%
% Burada:
%   mu = 398600 km³/s² = 3.986e14 m³/s²
%   r = km
%   I = kg·m²
%   Sonuç: N·m

% Birim dönüşümü
mu_m = mu * 1e9;        % km³/s² -> m³/s²
r_m = r_mag * 1e3;      % km -> m

% Atalet matrisi çarpımı
In = I * n_body;

% Çapraz çarpım
T_gg = (3 * mu_m / r_m^3) * cross(n_body, In);

end
