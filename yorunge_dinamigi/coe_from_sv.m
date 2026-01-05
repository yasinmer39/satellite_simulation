% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Appendix D.8
% Algorithm 4.1 sayfa 606

% Bu fonksiyon, state vektöründen (r, v) klasik yörüngesel elementleri hesaplar.

% mu   - kütleçekimsel parametre (km^3/s^2)
% R    - pozisyon vektörü (km)
% V    - hız vektörü (km/s)
% r    - pozisyon vektörünün büyüklüğü (km)
% v    - hız vektörünün büyüklüğü (km/s)
% vr   - radyal hız bileşeni (km/s)
% H    - açısal momentum vektörü (km^2/s)
% h    - açısal momentum vektörünün büyüklüğü (km^2/s)
% incl - eğilim, inclination (rad)
% N    - düğüm çizgisi vektörü (km^2/s)
% n    - düğüm çizgisi vektörünün büyüklüğü (km^2/s)
% RA   - yükselen düğümün sağ açıklığı, RAAN (rad)
% E    - dış merkezlilik vektörü, eccentricity vector
% e    - dış merkezlilik, eccentricity
% w    - yerberi açısı, argument of perigee (rad)
% TA   - gerçek anomali, true anomaly (rad)
% a    - yarı büyük eksen, semi-major axis (km)
% coe  - yörüngesel elementler vektörü [h e RA incl w TA a]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function coe = coe_from_sv(R, V)

global mu

eps = 1.e-10;  % küçük bir sayı tolerans için

r = norm(R);
v = norm(V);

% Radyal hız bileşeni (4.7 denklemi)
vr = dot(R, V) / r;

% Açısal momentum vektörü (4.8 denklemi)
H = cross(R, V);
h = norm(H);

% Eğilim (4.9 denklemi)
incl = acos(H(3) / h);

% Düğüm çizgisi vektörü (4.10 denklemi)
N = cross([0 0 1], H);
n = norm(N);

% Yükselen düğümün sağ açıklığı, RAAN (4.11 denklemi)
if n ~= 0
    RA = acos(N(1) / n);
    if N(2) < 0
        RA = 2*pi - RA;
    end
else
    RA = 0;
end

% Dış merkezlilik vektörü (4.12 denklemi)
E = 1/mu * ((v^2 - mu/r)*R - r*vr*V);
e = norm(E);

% Yerberi açısı (4.13 denklemi)
if n ~= 0
    if e > eps
        w = acos(dot(N, E) / (n*e));
        if E(3) < 0
            w = 2*pi - w;
        end
    else
        w = 0;
    end
else
    w = 0;
end

% Gerçek anomali (4.14 denklemi)
if e > eps
    TA = acos(dot(E, R) / (e*r));
    if vr < 0
        TA = 2*pi - TA;
    end
else
    cp = cross(N, R);
    if cp(3) >= 0
        TA = acos(dot(N, R) / (n*r));
    else
        TA = 2*pi - acos(dot(N, R) / (n*r));
    end
end

% Yarı büyük eksen (4.60 denklemi, eliptik yörüngeler için)
a = h^2 / mu / (1 - e^2);

% Yörüngesel elementler vektörü
coe = [h e RA incl w TA a];

end
