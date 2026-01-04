% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Appendix D.6
% Sayfa 603

% Bu fonksiyon, Lagrange f ve g katsayılarını hesaplar.

% x   - evrensel anomali x = 0'dan beri (km^0.5)
% t   - x = 0'dan beri geçen zaman (s)
% ro  - x = 0'da radyal pozisyon (km)
% a   - semi-major eksenin tersi (1/km)
% f   - f Lagrange katsayısı (boyutsuz)
% g   - g Lagrange katsayısı (s)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: stumpC, stumpS

function [f, g] = f_and_g(x, t, ro, a)

global mu

z = a*x^2;

% 3.66a denklemi
f = 1 - x^2/ro*stumpC(z);

% 3.66b denklemi
g = t - 1/sqrt(mu)*x^3*stumpS(z);

end
