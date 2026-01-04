% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Appendix D.4
% Sayfa 600

% Bu fonksiyon, Stumpff fonksiyonu S(z)'yi hesaplar.

% z - girdi argümanı
% s - S(z) Stumpff fonksiyonunun değeri

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function s = stumpS(z)

if z > 0
    s = (sqrt(z) - sin(sqrt(z))) / (sqrt(z))^3;
elseif z < 0
    s = (sinh(sqrt(-z)) - sqrt(-z)) / (sqrt(-z))^3;
else
    s = 1/6;
end

end
