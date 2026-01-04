% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Appendix D.4
% Sayfa 600

% Bu fonksiyon, Stumpff fonksiyonu C(z)'yi hesaplar.

% z - girdi argümanı
% c - C(z) Stumpff fonksiyonunun değeri

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function c = stumpC(z)

if z > 0
    c = (1 - cos(sqrt(z))) / z;
elseif z < 0
    c = (cosh(sqrt(-z)) - 1) / (-z);
else
    c = 1/2;
end

end
