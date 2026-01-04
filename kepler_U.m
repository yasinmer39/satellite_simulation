% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Appendix D.5
% Algorithm 3.3 sayfa 601

% Bu fonksiyon, evrensel anomali için evrensel Kepler denklemini çözmek için 
% Newton'un yöntemini kullanır.

% mu   - kütleçekimsel parametre (km^3/s^2)
% x    - evrensel anomali, universal anomaly (km^0.5)
% dt   - x = 0'dan beri geçen zaman (s)
% ro   - x = 0'da radyal pozisyon (km)
% vro  - x = 0'da radyal hız (km/s)
% a    - semi-major eksenin tersi (1/km)
% z    - ek değişken (z = a*x^2)
% C    - C(z) Stumpff fonksiyonunun değeri 
% S    - S(z) Stumpff fonksiyonunun değeri
% n    - yakınsama için iterasyon sayısı
% nMax - izin verilen maksimum iterasyon sayısı

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: stumpC, stumpS

function x = kepler_U(dt, ro, vro, a)

global mu

% Hata toleransı ve maksimum iterasyon sayısını ayarla
error_tol = 1.e-8;
nMax = 1000;

% x için başlangıç tahmini (3.48 denklemi)
x = sqrt(mu) * abs(a) * dt;

% 3.45-3.47 denklemlerini kullanarak Newton iterasyonu
n = 0;
ratio = 1;

while abs(ratio) > error_tol && n <= nMax
    n = n + 1;
    C = stumpC(a*x^2);
    S = stumpS(a*x^2);
    F = ro*vro/sqrt(mu)*x^2*C + (1 - a*ro)*x^3*S + ro*x - sqrt(mu)*dt;
    dFdx = ro*vro/sqrt(mu)*x*(1 - a*x^2*S) + (1 - a*ro)*x^2*C + ro;
    ratio = F/dFdx;
    x = x - ratio;
end

% Yakınsama kontrolü
if n > nMax
    fprintf('\n **kepler_U %g iterasyonda yakınsamadı\n\n', nMax)
end

end
