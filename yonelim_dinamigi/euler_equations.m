% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Denklem 16-50: Euler's Dynamic Equations of Motion

% Bu fonksiyon, rijit cisim için Euler dinamik denklemlerini hesaplar.
% Principal (asal) eksenler koordinat sisteminde çalışır.

% omega     - açısal hız vektörü body frame'de (rad/s) [3x1]
% I         - asal atalet momentleri [I1, I2, I3] (kg·m²) [3x1] veya [1x3]
% N         - dış tork vektörü body frame'de (N·m) [3x1]
% omega_dot - açısal ivme vektörü body frame'de (rad/s²) [3x1]

% Euler Denklemleri (Wertz Eq. 16-50):
%   I1 * dω1/dt = N1 + (I2 - I3) * ω2 * ω3
%   I2 * dω2/dt = N2 + (I3 - I1) * ω3 * ω1
%   I3 * dω3/dt = N3 + (I1 - I2) * ω1 * ω2

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function omega_dot = euler_equations(omega, I, N)

% Varsayılan tork (sıfır)
if nargin < 3
    N = [0; 0; 0];
end

% Vektörleri sütun vektörüne çevir
if isrow(omega)
    omega = omega';
end
if isrow(I)
    I = I';
end
if isrow(N)
    N = N';
end

% Asal atalet momentleri
I1 = I(1);
I2 = I(2);
I3 = I(3);

% Açısal hız bileşenleri
w1 = omega(1);
w2 = omega(2);
w3 = omega(3);

% Euler denklemleri (Wertz Eq. 16-50)
omega_dot = zeros(3, 1);

omega_dot(1) = (N(1) + (I2 - I3) * w2 * w3) / I1;
omega_dot(2) = (N(2) + (I3 - I1) * w3 * w1) / I2;
omega_dot(3) = (N(3) + (I1 - I2) * w1 * w2) / I3;

end
