% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Appendix E - Coordinate Transformations

% Bu fonksiyon, 3-2-1 (Yaw-Pitch-Roll) Euler açılarından DCM hesaplar.

% euler     - Euler açıları [phi, theta, psi]' (rad) [3x1]
%             phi   = roll  (x ekseni etrafında) - 1. dönüş body'de
%             theta = pitch (y ekseni etrafında) - 2. dönüş  
%             psi   = yaw   (z ekseni etrafında) - 3. dönüş inertial'da
% A         - DCM [3x3] - Body'den Reference'a dönüşüm

% Dönüş sırası: R = R_z(psi) * R_y(theta) * R_x(phi)
%               (İnertial'dan body'ye: sıra tersine)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function A = euler2dcm_321(euler)

% Vektörü sütun vektörüne çevir
if isrow(euler)
    euler = euler';
end

% Euler açıları
phi = euler(1);     % roll
theta = euler(2);   % pitch
psi = euler(3);     % yaw

% Trigonometrik terimler
c_phi = cos(phi);     s_phi = sin(phi);
c_theta = cos(theta); s_theta = sin(theta);
c_psi = cos(psi);     s_psi = sin(psi);

% Temel dönüşüm matrisleri
R_x = [1    0       0    ;
       0    c_phi   s_phi;
       0   -s_phi   c_phi];

R_y = [c_theta   0  -s_theta;
       0         1   0      ;
       s_theta   0   c_theta];

R_z = [ c_psi   s_psi   0;
       -s_psi   c_psi   0;
        0       0       1];

% DCM: Body'den Reference'a
% A = (R_z * R_y * R_x)' = R_x' * R_y' * R_z'
% Ancak convention'a göre A = R_z * R_y * R_x (body'den ref'e)

A = R_z * R_y * R_x;

end
