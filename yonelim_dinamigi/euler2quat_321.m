% Spacecraft Attitude Determination and Control (Wertz) kitabından uyarlanmıştır.
% Section 12.1 - Parameterization of the Attitude

% Bu fonksiyon, 3-2-1 Euler açılarından quaternion hesaplar.

% euler     - Euler açıları [phi, theta, psi]' (rad) [3x1]
%             phi   = roll  (x ekseni etrafında)
%             theta = pitch (y ekseni etrafında)  
%             psi   = yaw   (z ekseni etrafında)
% q         - quaternion [q1, q2, q3, q4]' (q4 skaler kısım) [4x1]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function q = euler2quat_321(euler)

% Vektörü sütun vektörüne çevir
if isrow(euler)
    euler = euler';
end

% Euler açıları
phi = euler(1);     % roll
theta = euler(2);   % pitch
psi = euler(3);     % yaw

% Yarım açılar
c_phi2 = cos(phi/2);    s_phi2 = sin(phi/2);
c_theta2 = cos(theta/2); s_theta2 = sin(theta/2);
c_psi2 = cos(psi/2);    s_psi2 = sin(psi/2);

% Quaternion hesaplama
q = zeros(4, 1);

q(1) = s_phi2*c_theta2*c_psi2 - c_phi2*s_theta2*s_psi2;
q(2) = c_phi2*s_theta2*c_psi2 + s_phi2*c_theta2*s_psi2;
q(3) = c_phi2*c_theta2*s_psi2 - s_phi2*s_theta2*c_psi2;
q(4) = c_phi2*c_theta2*c_psi2 + s_phi2*s_theta2*s_psi2;

% Normalize et
q = q / norm(q);

end
