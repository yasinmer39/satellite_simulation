% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Appendix D.7
% Algorithm 3.4 sayfa 604

% Bu fonksiyon, initial state vektöründen (R0, V0) ve aradan geçen süreden, 
% state vektörünü (R, V) hesaplar.

% mu - kütleçekimsel parametre (km^3/s^2)
% R0 - initial pozisyon vektörü (km)
% V0 - initial hız vektörü (km/s)
% t  - aradan geçen süre (s)
% R  - hesaplanan pozisyon vektörü (km)
% V  - hesaplanan hız vektörü (km/s)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: kepler_U,
% f_and_g, fDot_and_gDot

function [R, V] = rv_from_r0v0(R0, V0, t)

global mu

% R0 ve V0'ın büyüklükleri
r0 = norm(R0);
v0 = norm(V0);

% Radyal hız bileşeni (3.21 denklemi)
vr0 = dot(R0, V0) / r0;

% Semi-major eksenin tersi (3.44 denklemi)
alpha = 2/r0 - v0^2/mu;

% Evrensel anomaliyi hesapla (Algorithm 3.3)
x = kepler_U(t, r0, vr0, alpha);

% f ve g Lagrange katsayılarını hesapla (3.66a ve 3.66b denklemleri)
[f, g] = f_and_g(x, t, r0, alpha);

% t anındaki pozisyon vektörünü hesapla (3.64 denklemi)
R = f*R0 + g*V0;

% t anındaki pozisyonun büyüklüğü
r = norm(R);

% fdot ve gdot Lagrange katsayılarını hesapla (3.66c ve 3.66d denklemleri)
[fdot, gdot] = fDot_and_gDot(x, r, r0, alpha);

% t anındaki hız vektörünü hesapla (3.65 denklemi)
V = fdot*R0 + gdot*V0;

end
