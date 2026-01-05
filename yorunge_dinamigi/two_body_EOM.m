% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Chapter 2 - The Two-Body Problem
% Denklem 2.15 (sayfa 36)

% Bu fonksiyon, iki cisim probleminin hareket denklemini tanımlar.
% ODE çözücüler (ode45, ode113 vb.) ile kullanılmak üzere tasarlanmıştır.

% t     - zaman (s) - ODE çözücü tarafından sağlanır
% state - state vektörü [x; y; z; vx; vy; vz] (km ve km/s)
% mu    - kütleçekimsel parametre (km^3/s^2)
% dstate - state vektörünün türevi [vx; vy; vz; ax; ay; az]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function dstate = two_body_EOM(t, state)

global mu

% State vektöründen pozisyon ve hız çıkar
r = state(1:3);     % pozisyon vektörü (km)
v = state(4:6);     % hız vektörü (km/s)

% Pozisyon vektörünün büyüklüğü
r_mag = norm(r);

% İki cisim ivmesi (merkezi çekim kuvveti)
a = -mu / r_mag^3 * r;

% State türevi: [hız; ivme]
dstate = [v; a];

end
