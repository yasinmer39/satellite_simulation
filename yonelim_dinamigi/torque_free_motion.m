% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Denklem 16-59 ve 16-60: Torque-Free Motion (Axisymmetric Case)

% Bu fonksiyon, eksenel simetrik (I1 = I2) rijit cisim için
% torksuz hareketi analitik olarak hesaplar.

% t         - zaman (s) [skaler veya vektör]
% omega0    - başlangıç açısal hızı body frame'de (rad/s) [3x1]
% I         - asal atalet momentleri [I1, I2, I3] (kg·m²) [3x1]
%             I1 = I2 (eksenel simetri gerekli)
% omega     - açısal hız zaman serileri [3 x length(t)]
% omega_p   - body nutation rate (rad/s)

% Torque-Free Motion (Wertz Eq. 16-60):
%   ω1(t) = ω01*cos(ωp*t) + ω02*sin(ωp*t)
%   ω2(t) = ω02*cos(ωp*t) - ω01*sin(ωp*t)
%   ω3(t) = ω03 (sabit)
%
% Body nutation rate (Wertz Eq. 16-59b):
%   ωp = (I3 - IT)/IT * ω3

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function [omega, omega_p] = torque_free_motion(t, omega0, I)

% Vektörleri sütun vektörüne çevir
if isrow(omega0)
    omega0 = omega0';
end
if isrow(I)
    I = I';
end

% Asal atalet momentleri
I_T = I(1);     % Transverse moment (I1 = I2)
I_3 = I(3);     % Symmetry axis moment

% Eksenel simetri kontrolü
if abs(I(1) - I(2)) > 0.01 * I(1)
    warning('I1 ≠ I2. Bu fonksiyon eksenel simetri varsayar.');
end

% Başlangıç açısal hız bileşenleri
w01 = omega0(1);
w02 = omega0(2);
w03 = omega0(3);

% Body nutation rate (Wertz Eq. 16-59b)
omega_p = (I_3 - I_T) / I_T * w03;

% Zaman vektörünü satır vektörüne çevir
if iscolumn(t)
    t = t';
end

% Açısal hız zaman serileri (Wertz Eq. 16-60)
omega = zeros(3, length(t));

omega(1,:) = w01 * cos(omega_p * t) + w02 * sin(omega_p * t);
omega(2,:) = w02 * cos(omega_p * t) - w01 * sin(omega_p * t);
omega(3,:) = w03 * ones(1, length(t));

end
