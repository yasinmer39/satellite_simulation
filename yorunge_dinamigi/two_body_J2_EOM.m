% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Chapter 4.7 - Effects of Earth's oblateness
% İki cisim + J2 pertürbasyonu

% Bu fonksiyon, J2 pertürbasyonu dahil hareket denklemini tanımlar.
% ODE çözücüler (ode45, ode113 vb.) ile kullanılmak üzere tasarlanmıştır.

% t     - zaman (s) - ODE çözücü tarafından sağlanır
% state - state vektörü [x; y; z; vx; vy; vz] (km ve km/s)
% mu    - kütleçekimsel parametre (km^3/s^2)
% dstate - state vektörünün türevi [vx; vy; vz; ax; ay; az]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: J2_perturbation

function dstate = two_body_J2_EOM(t, state)

global mu

% State vektöründen pozisyon ve hız çıkar
r = state(1:3);     % pozisyon vektörü (km)
v = state(4:6);     % hız vektörü (km/s)

% Pozisyon vektörünün büyüklüğü
r_mag = norm(r);

% İki cisim ivmesi (merkezi çekim kuvveti)
a_two_body = -mu / r_mag^3 * r;

% J2 pertürbasyon ivmesi
a_J2 = J2_perturbation(r);

% Toplam ivme
a_total = a_two_body + a_J2;

% State türevi: [hız; ivme]
dstate = [v; a_total];

end
