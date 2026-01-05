% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Denklem 16-34: Angular Momentum

% Bu fonksiyon, açısal momentum vektörünü hesaplar.

% omega     - açısal hız vektörü body frame'de (rad/s) [3x1]
% I         - atalet matrisi veya asal momentler [3x3] veya [3x1]
% L         - açısal momentum vektörü body frame'de (kg·m²/s) [3x1]

% Angular Momentum (Wertz Eq. 16-34):
%   L = I * ω

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function L = angular_momentum(omega, I)

% Vektörü sütun vektörüne çevir
if isrow(omega)
    omega = omega';
end

% I matris mi vektör mü kontrol et
if isvector(I)
    % Asal momentler verilmiş, diagonal matris oluştur
    if isrow(I)
        I = I';
    end
    I_matrix = diag(I);
else
    I_matrix = I;
end

% Angular momentum
L = I_matrix * omega;

end
