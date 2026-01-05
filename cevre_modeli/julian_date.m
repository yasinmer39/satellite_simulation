% Orbital Mechanics for Engineering Students (Curtis) kitabından alınmıştır.
% Algorithm 5.2 - Julian Date

% Bu fonksiyon, Gregoryen takvim tarihinden Julian Date hesaplar.

% year   - yıl (örn: 2026)
% month  - ay (1-12)
% day    - gün (1-31)
% hour   - saat (0-23)
% minute - dakika (0-59)
% second - saniye (0-59)
% JD     - Julian Date

% Örnek: 1 Ocak 2000, 12:00 UT -> JD = 2451545.0

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function JD = julian_date(year, month, day, hour, minute, second)

% Varsayılan değerler
if nargin < 4
    hour = 0;
end
if nargin < 5
    minute = 0;
end
if nargin < 6
    second = 0;
end

% Gün kesri
UT = hour + minute/60 + second/3600;    % Saat cinsinden

% Julian Date algoritması (Curtis Denklem 5.46)
J0 = 367*year - floor(7*(year + floor((month + 9)/12))/4) ...
     + floor(275*month/9) + day + 1721013.5;

JD = J0 + UT/24;

end
