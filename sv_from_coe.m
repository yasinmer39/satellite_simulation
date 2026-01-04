% Orbital Mechanics for Engineering Students kitabından alınmıştır.
% Appendix D.9
% Algorithm 4.2 sayfa 610

% Bu fonksiyon, (r,v) state vektörünü, klasik yörüngesel elementlerini kullanarak hesaplar.

% mu   - kütleçekimsel parametre (km^3/s^2)
% coe  - yörüngesel elementler [h e RA incl w TA]:
%       h    = açısal momentum (km^2/s)
%       e    = dış merkezlilik, eccentricity
%       RA   = yükselen düğümün sağ açıklığı, RAAN (rad)
%       incl = eğilim, inclination of the orbit (rad)
%       w    = yerberi açısı, argument of perigee (rad)
%       TA   = gerçek anomali, true anomaly (rad)
% R3_w - z-ekseni etrafında w açısı ile yapılan rotasyon matrisi
% R1_i - x-ekseni etrafında i açısı ile yapılan rotasyon matrisi
% R3_W - z-ekseni etrafında RA açısı ile yapılan rotasyon matrisi
% Q_pX - perifokal koordinat sisteminden yermerkezli ekvator koordinat sistemine dönüşüm matrisi
% rp   - perifokal koordinat sistemindeki konum vektörü (km)
% vp   - perifokal koordinat sistemindeki hız vektörü (km/s)
% r    - yermerkezli ekvator koordinat sistemindeki konum vektörü (km)
% v    - yermerkezli ekvator koordinat sistemindeki hız vektörü (km/s)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [r, v] = sv_from_coe(coe)

global mu

h    = coe(1);
e    = coe(2);
RA   = coe(3);
incl = coe(4);
w    = coe(5);
TA   = coe(6);

% Perifokal koordinat sisteminde pozisyon ve hız
% Kitaptaki 4.37 ve 4.38 denklemleri
rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

% Rotasyon matrisleri
% 4.39 denklemi - z ekseni etrafında RA (RAAN) açısı ile rotasyon
R3_W = [ cos(RA) sin(RA) 0
        -sin(RA) cos(RA) 0
            0       0    1];

% 4.40 denklemi - x ekseni etrafında incl (eğilim) açısı ile rotasyon
R1_i = [ 1    0          0
         0 cos(incl)  sin(incl)
         0 -sin(incl) cos(incl)];

% 4.41 denklemi - z ekseni etrafında w (argument of perigee) açısı ile rotasyon
R3_w = [ cos(w) sin(w)   0
        -sin(w) cos(w)   0
            0       0    1];

% 4.44 denklemi - Perifokal'dan Geocentric Equatorial'a dönüşüm matrisi
Q_pX = R3_W' * R1_i' * R3_w';

% 4.46 denklemi - State vektörlerini dönüştür (r ve v sütun vektörleridir)
r = Q_pX * rp;
v = Q_pX * vp;

% r ve v vektörlerini satır vektörlerine çevir
r = r';
v = v';

end
