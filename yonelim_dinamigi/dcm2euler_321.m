% Spacecraft Attitude Determination and Control (Wertz) kitabından uyarlanmıştır.
% Appendix E - Coordinate Transformations

% Bu fonksiyon, DCM'den 3-2-1 (Yaw-Pitch-Roll) Euler açıları hesaplar.

% A         - DCM [3x3]
% euler     - Euler açıları [phi, theta, psi]' (rad) [3x1]
%             phi   = roll  (x ekseni etrafında)
%             theta = pitch (y ekseni etrafında)  
%             psi   = yaw   (z ekseni etrafında)

% Singularity: theta = ±90° (gimbal lock)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function euler = dcm2euler_321(A)

euler = zeros(3, 1);

% Pitch (theta) - asin'den
% A(1,3) = -sin(theta)
sin_theta = -A(1,3);

% Singularity kontrolü
if abs(sin_theta) > 0.99999
    % Gimbal lock durumu
    warning('Gimbal lock! theta ≈ ±90°');
    
    % theta = ±90°
    euler(2) = sign(sin_theta) * pi/2;
    
    % Bu durumda phi ve psi bağımlı
    % Sadece (phi - psi) veya (phi + psi) belirlenebilir
    euler(1) = 0;  % phi = 0 varsayımı
    euler(3) = atan2(-A(2,1), A(2,2));  % psi
else
    % Normal durum
    euler(2) = asin(sin_theta);  % theta
    
    % Roll (phi)
    % A(2,3) = sin(phi) * cos(theta)
    % A(3,3) = cos(phi) * cos(theta)
    euler(1) = atan2(A(2,3), A(3,3));  % phi
    
    % Yaw (psi)
    % A(1,2) = sin(psi) * cos(theta)
    % A(1,1) = cos(psi) * cos(theta)
    euler(3) = atan2(A(1,2), A(1,1));  % psi
end

end
