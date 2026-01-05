% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Section 5.5 - Modeling the Positions of the Sun, Moon, and Planets
% Orbital Mechanics for Engineering Students (Curtis) Algorithm 5.3 ile desteklenmiştir.

% Bu fonksiyon, Güneş'in ECI (Geocentric Equatorial Inertial) koordinatlarındaki
% pozisyonunu hesaplar.

% JD       - Julian Date
% r_sun    - Güneş pozisyon vektörü ECI frame'de (km) [3x1]
% lambda   - Güneş ekliptik boylamı (rad)
% epsilon  - Ekliptik eğikliği (rad)

% Doğruluk: ~0.01 derece (attitude determination için yeterli)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [r_sun, lambda, epsilon] = sun_position(JD)

% Güneş-Dünya ortalama mesafesi (1 AU)
AU = 149597870.7;   % km

% J2000.0 epoch'tan beri geçen Julian yüzyıllar
% J2000.0 = JD 2451545.0 (1 Ocak 2000, 12:00 TT)
T = (JD - 2451545.0) / 36525;

% Güneş'in ortalama boylamı (derece)
% Wertz Denklem 5-43 benzeri
L0 = 280.4606184 + 36000.77005361 * T;
L0 = mod(L0, 360);
if L0 < 0
    L0 = L0 + 360;
end

% Güneş'in ortalama anomalisi (derece)
M = 357.5291092 + 35999.05034 * T;
M = mod(M, 360);
if M < 0
    M = M + 360;
end

% Dereceyi radyana çevir
M_rad = M * pi / 180;

% Ekliptik boylam (derece) - center denklemi ile düzeltilmiş
% Curtis Denklem 5.48
lambda_deg = L0 + 1.914666471 * sin(M_rad) + 0.019994643 * sin(2*M_rad);
lambda = lambda_deg * pi / 180;     % radyan

% Ekliptik eğikliği (obliquity of the ecliptic) - derece
% Curtis Denklem 5.49
epsilon_deg = 23.439291 - 0.0130042 * T;
epsilon = epsilon_deg * pi / 180;   % radyan

% Güneş mesafesi (AU)
% Curtis Denklem 5.50
r_AU = 1.000140612 - 0.016708617 * cos(M_rad) - 0.000139589 * cos(2*M_rad);
r_km = r_AU * AU;                   % km

% ECI koordinatları
% Güneş ekliptik düzlemde, sonra equatorial'e dönüştür
% r_sun = R1(-epsilon) * [r*cos(lambda); r*sin(lambda); 0]
cos_lambda = cos(lambda);
sin_lambda = sin(lambda);
cos_epsilon = cos(epsilon);
sin_epsilon = sin(epsilon);

r_sun = r_km * [cos_lambda;
                sin_lambda * cos_epsilon;
                sin_lambda * sin_epsilon];

end
