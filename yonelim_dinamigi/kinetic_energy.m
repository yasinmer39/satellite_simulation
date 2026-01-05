% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Denklem 16-39: Rotational Kinetic Energy

% Bu fonksiyon, dönme kinetik enerjisini hesaplar.

% omega     - açısal hız vektörü body frame'de (rad/s) [3x1]
% I         - atalet matrisi veya asal momentler [3x3] veya [3x1]
% Ek        - dönme kinetik enerjisi (J = kg·m²/s²)

% Kinetik Enerji (Wertz Eq. 16-39):
%   Ek = (1/2) * ω' * I * ω
%      = (1/2) * (I1*ω1² + I2*ω2² + I3*ω3²)  [asal eksenlerde]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function Ek = kinetic_energy(omega, I)

% Vektörü sütun vektörüne çevir
if isrow(omega)
    omega = omega';
end

% I matris mi vektör mü kontrol et
if isvector(I)
    % Asal momentler verilmiş
    if isrow(I)
        I = I';
    end
    Ek = 0.5 * (I(1)*omega(1)^2 + I(2)*omega(2)^2 + I(3)*omega(3)^2);
else
    % Tam atalet matrisi
    Ek = 0.5 * omega' * I * omega;
end

end
