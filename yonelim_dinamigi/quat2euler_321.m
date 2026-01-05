% Spacecraft Attitude Determination and Control (Wertz) kitabından uyarlanmıştır.
% Section 12.1 - Parameterization of the Attitude

% Bu fonksiyon, quaternion'dan 3-2-1 Euler açıları hesaplar.

% q         - quaternion [q1, q2, q3, q4]' (q4 skaler kısım) [4x1]
% euler     - Euler açıları [phi, theta, psi]' (rad) [3x1]
%             phi   = roll  (x ekseni etrafında)
%             theta = pitch (y ekseni etrafında)  
%             psi   = yaw   (z ekseni etrafında)

% Singularity: theta = ±90° (gimbal lock)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function euler = quat2euler_321(q)

% Vektörü sütun vektörüne çevir
if isrow(q)
    q = q';
end

% Normalize et
q = q / norm(q);

% Quaternion bileşenleri
q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);

euler = zeros(3, 1);

% DCM elemanlarını quaternion cinsinden hesapla
% A(1,3) = 2*(q1*q3 - q2*q4) = -sin(theta)
% A(2,3) = 2*(q2*q3 + q1*q4) = sin(phi)*cos(theta)
% A(3,3) = q4^2 - q1^2 - q2^2 + q3^2 = cos(phi)*cos(theta)
% A(1,2) = 2*(q1*q2 + q3*q4) = sin(psi)*cos(theta)
% A(1,1) = q4^2 + q1^2 - q2^2 - q3^2 = cos(psi)*cos(theta)

% Pitch (theta)
sin_theta = 2 * (q2*q4 - q1*q3);

% Singularity kontrolü
if abs(sin_theta) > 0.99999
    % Gimbal lock
    warning('Gimbal lock! theta ≈ ±90°');
    
    euler(2) = sign(sin_theta) * pi/2;
    euler(1) = 0;
    euler(3) = 2 * atan2(q1, q4);
else
    % Normal durum
    euler(2) = asin(sin_theta);  % theta
    
    % Roll (phi)
    euler(1) = atan2(2*(q3*q4 + q1*q2), ...
                     q4^2 - q1^2 - q2^2 + q3^2);
    
    % Yaw (psi)
    euler(3) = atan2(2*(q1*q4 + q2*q3), ...
                     q4^2 + q1^2 - q2^2 - q3^2);
end

end
