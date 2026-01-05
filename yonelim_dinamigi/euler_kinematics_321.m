% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 ve Appendix E - Euler Angle Kinematics
% 3-2-1 (Yaw-Pitch-Roll) dizisi için kinematik denklemler

% Bu fonksiyon, 3-2-1 Euler açıları için kinematik denklemleri hesaplar.
% DİKKAT: theta = ±90° civarında singularity (gimbal lock) vardır!

% euler     - Euler açıları [phi, theta, psi]' (rad) [3x1]
%             phi   = roll  (x ekseni etrafında)
%             theta = pitch (y ekseni etrafında)  
%             psi   = yaw   (z ekseni etrafında)
% omega     - açısal hız vektörü body frame'de (rad/s) [3x1]
% euler_dot - Euler açı türevleri (rad/s) [3x1]

% 3-2-1 Euler Kinematik:
%   φ̇ = ω1 + (ω2 sin(φ) + ω3 cos(φ)) tan(θ)
%   θ̇ = ω2 cos(φ) - ω3 sin(φ)
%   ψ̇ = (ω2 sin(φ) + ω3 cos(φ)) / cos(θ)

% Singularity: θ = ±90° (cos(θ) = 0)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function euler_dot = euler_kinematics_321(euler, omega)

% Vektörleri sütun vektörüne çevir
if isrow(euler)
    euler = euler';
end
if isrow(omega)
    omega = omega';
end

% Euler açıları
phi = euler(1);     % roll
theta = euler(2);   % pitch
psi = euler(3);     % yaw (kullanılmıyor ama tutarlılık için)

% Açısal hız bileşenleri
w1 = omega(1);  % p (roll rate)
w2 = omega(2);  % q (pitch rate)
w3 = omega(3);  % r (yaw rate)

% Trigonometrik terimler
sin_phi = sin(phi);
cos_phi = cos(phi);
tan_theta = tan(theta);
cos_theta = cos(theta);

% Singularity kontrolü
if abs(cos_theta) < 1e-10
    warning('Gimbal lock! theta ≈ ±90°. Quaternion kullanın.');
    cos_theta = sign(cos_theta) * 1e-10;
end

% Euler kinematik denklemleri (3-2-1 dizisi)
euler_dot = zeros(3, 1);

euler_dot(1) = w1 + (w2 * sin_phi + w3 * cos_phi) * tan_theta;  % φ̇
euler_dot(2) = w2 * cos_phi - w3 * sin_phi;                      % θ̇
euler_dot(3) = (w2 * sin_phi + w3 * cos_phi) / cos_theta;        % ψ̇

end
