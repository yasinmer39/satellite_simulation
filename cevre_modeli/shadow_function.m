% Spacecraft Attitude Determination and Control (Wertz) kitabından uyarlanmıştır.
% Section 5.3 ve Orbital Mechanics for Engineering Students (Curtis) Section 10.7

% Bu fonksiyon, uydunun Dünya'nın gölgesinde olup olmadığını hesaplar.
% Hem umbra (tam gölge) hem de penumbra (yarı gölge) durumlarını kontrol eder.

% r_sat     - uydu pozisyonu ECI frame'de (km) [3x1]
% r_sun     - güneş pozisyonu ECI frame'de (km) [3x1]
% nu        - gölge faktörü (0 = tam gölge, 1 = tam güneş, 0-1 arası = penumbra)
% in_shadow - gölgede mi? (boolean, umbra için)

% Geometri:
%   - Dünya yarıçapı: R_earth
%   - Güneş yarıçapı: R_sun
%   - Umbra: Güneş tamamen Dünya arkasında
%   - Penumbra: Güneş kısmen Dünya arkasında

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [nu, in_shadow] = shadow_function(r_sat, r_sun)

% Sabitler
R_earth = 6378.137;         % Dünya yarıçapı (km)
R_sun = 696000;             % Güneş yarıçapı (km)

% Pozisyon vektörlerini sütun vektörüne çevir
if isrow(r_sat)
    r_sat = r_sat';
end
if isrow(r_sun)
    r_sun = r_sun';
end

% Vektör büyüklükleri
r_sat_mag = norm(r_sat);
r_sun_mag = norm(r_sun);

% Güneş yön birim vektörü (Dünya merkezinden)
s_hat = r_sun / r_sun_mag;

% Uydu pozisyonunun güneş yönündeki projeksiyonu
% r_sat · s_hat = r_sat_mag * cos(theta)
r_dot_s = dot(r_sat, s_hat);

% Gölge kontrolü için geometri
% Uydu güneşin arkasındaysa (r_dot_s > 0), gölge mümkün değil
if r_dot_s > 0
    nu = 1.0;
    in_shadow = false;
    return;
end

% Uydunun güneş yönüne dik mesafesi
% d = |r_sat - (r_sat · s_hat) * s_hat|
r_perp = r_sat - r_dot_s * s_hat;
d = norm(r_perp);

% Umbra yarıçapı (uydunun güneşe olan mesafesinde)
% Basitleştirilmiş silindirik gölge modeli
% Gerçek konik model için daha karmaşık hesaplama gerekir
r_umbra = R_earth;          % Basit silindirik model

% Penumbra yarıçapı
% Konik modelde: r_penumbra = R_earth + (R_sun + R_earth) * |r_dot_s| / r_sun_mag
r_penumbra = R_earth * 1.02;    % Basit yaklaşım (~2% daha geniş)

% Gölge durumu kontrolü
if d < r_umbra
    % Umbra (tam gölge)
    nu = 0.0;
    in_shadow = true;
elseif d < r_penumbra
    % Penumbra (yarı gölge)
    % Lineer interpolasyon
    nu = (d - r_umbra) / (r_penumbra - r_umbra);
    in_shadow = false;      % Teknik olarak penumbra
else
    % Güneşte
    nu = 1.0;
    in_shadow = false;
end

end
