% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Denklem 16-54: Spacecraft with Reaction Wheels

% Bu fonksiyon, reaksiyon çarklı uydu için Euler denklemlerini hesaplar.

% omega     - uydu açısal hızı body frame'de (rad/s) [3x1]
% I         - uydu atalet matrisi (çarklar dahil) [3x3] veya [3x1]
% h_w       - çark açısal momentum vektörü body frame'de (N·m·s) [3x1]
% h_w_dot   - çark açısal momentum türevi (uygulanan çark torku) [3x1]
% N_ext     - dış tork vektörü (N·m) [3x1]
% omega_dot - uydu açısal ivmesi (rad/s²) [3x1]

% Euler Denklemleri (Reaction Wheel ile) - Wertz Eq. 16-54:
%   I * dω/dt = N_ext - dh_w/dt - ω × (I*ω + h_w)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function omega_dot = euler_equations_with_wheels(omega, I, h_w, h_w_dot, N_ext)

% Varsayılan değerler
if nargin < 4
    h_w_dot = [0; 0; 0];
end
if nargin < 5
    N_ext = [0; 0; 0];
end

% Vektörleri sütun vektörüne çevir
if isrow(omega)
    omega = omega';
end
if isrow(h_w)
    h_w = h_w';
end
if isrow(h_w_dot)
    h_w_dot = h_w_dot';
end
if isrow(N_ext)
    N_ext = N_ext';
end

% I matris mi vektör mü kontrol et
if isvector(I)
    if isrow(I)
        I = I';
    end
    I_matrix = diag(I);
else
    I_matrix = I;
end

% Toplam açısal momentum (uydu + çarklar)
L_total = I_matrix * omega + h_w;

% Euler denklemi (Wertz Eq. 16-54)
% I * dω/dt = N_ext - dh_w/dt - ω × L_total
rhs = N_ext - h_w_dot - cross(omega, L_total);

% Açısal ivme
omega_dot = I_matrix \ rhs;

end
