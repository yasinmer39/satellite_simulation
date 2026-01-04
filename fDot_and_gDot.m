% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Appendix D.6
% Sayfa 603

% Bu fonksiyon, Lagrange f ve g katsayılarının zaman türevlerini hesaplar.

% x    - evrensel anomali x = 0'dan beri (km^0.5)
% r    - t anında radyal pozisyon (km)
% ro   - x = 0'da radyal pozisyon (km)
% a    - semi-major eksenin tersi (1/km)
% fdot - f'nin zaman türevi (1/s)
% gdot - g'nin zaman türevi (boyutsuz)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: stumpC, stumpS

function [fdot, gdot] = fDot_and_gDot(x, r, ro, a)

global mu

z = a*x^2;

% 3.66c denklemi
fdot = sqrt(mu)/(r*ro) * (a*x^3*stumpS(z) - x);

% 3.66d denklemi
gdot = 1 - x^2/r*stumpC(z);

end
