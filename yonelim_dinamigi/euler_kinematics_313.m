% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Denklem 16-12: 3-1-3 Euler Angle Kinematics

% Bu fonksiyon, 3-1-3 Euler açıları için kinematik denklemleri hesaplar.
% DİKKAT: theta = 0° veya 180° civarında singularity vardır!

% euler     - Euler açıları [phi, theta, psi]' (rad) [3x1]
%             phi   = ilk z ekseni etrafında dönüş
%             theta = x' ekseni etrafında dönüş
%             psi   = son z'' ekseni etrafında dönüş
% omega     - açısal hız vektörü body frame'de (rad/s) [3x1]
% euler_dot - Euler açı türevleri (rad/s) [3x1]

% 3-1-3 Euler Kinematik (Wertz Eq. 16-12):
%   θ̇ = ω_u cos(ψ) - ω_v sin(ψ)
%   φ̇ = (ω_u sin(ψ) + ω_v cos(ψ)) / sin(θ)
%   ψ̇ = ω_w - (ω_u sin(ψ) + ω_v cos(ψ)) cot(θ)

% Singularity: θ = 0° veya 180° (sin(θ) = 0)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function euler_dot = euler_kinematics_313(euler, omega)

% Vektörleri sütun vektörüne çevir
if isrow(euler)
    euler = euler';
end
if isrow(omega)
    omega = omega';
end

% Euler açıları
phi = euler(1);     % φ
theta = euler(2);   % θ
psi = euler(3);     % ψ

% Açısal hız bileşenleri (body frame: u, v, w)
w_u = omega(1);
w_v = omega(2);
w_w = omega(3);

% Trigonometrik terimler
sin_psi = sin(psi);
cos_psi = cos(psi);
sin_theta = sin(theta);
cot_theta = cot(theta);

% Singularity kontrolü
if abs(sin_theta) < 1e-10
    warning('Gimbal lock! theta ≈ 0° veya 180°. Quaternion kullanın.');
    sin_theta = sign(sin_theta) * 1e-10;
    cot_theta = cos(theta) / sin_theta;
end

% Euler kinematik denklemleri (Wertz Eq. 16-12)
euler_dot = zeros(3, 1);

euler_dot(1) = (w_u * sin_psi + w_v * cos_psi) / sin_theta;                 % φ̇
euler_dot(2) = w_u * cos_psi - w_v * sin_psi;                                % θ̇
euler_dot(3) = w_w - (w_u * sin_psi + w_v * cos_psi) * cot_theta;           % ψ̇

end
