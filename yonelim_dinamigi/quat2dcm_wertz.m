% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Section 12.1 - Parameterization of the Attitude
% Denklem 12-14: Quaternion to DCM Conversion

% Bu fonksiyon, quaternion'dan Direction Cosine Matrix (DCM) hesaplar.

% q         - quaternion [q1, q2, q3, q4]' (q4 skaler kısım) [4x1]
% A         - DCM (yön kosinüsleri matrisi) [3x3]
%             Body frame'den Reference frame'e dönüşüm

% DCM formülü (Wertz Eq. 12-14):
%   A = [q4² + q1² - q2² - q3²    2(q1q2 + q3q4)         2(q1q3 - q2q4)     ]
%       [2(q1q2 - q3q4)           q4² - q1² + q2² - q3²   2(q2q3 + q1q4)     ]
%       [2(q1q3 + q2q4)           2(q2q3 - q1q4)          q4² - q1² - q2² + q3²]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function A = quat2dcm_wertz(q)

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
q4 = q(4);  % skaler kısım

% DCM hesaplama (Wertz Eq. 12-14)
A = zeros(3, 3);

A(1,1) = q4^2 + q1^2 - q2^2 - q3^2;
A(1,2) = 2 * (q1*q2 + q3*q4);
A(1,3) = 2 * (q1*q3 - q2*q4);

A(2,1) = 2 * (q1*q2 - q3*q4);
A(2,2) = q4^2 - q1^2 + q2^2 - q3^2;
A(2,3) = 2 * (q2*q3 + q1*q4);

A(3,1) = 2 * (q1*q3 + q2*q4);
A(3,2) = 2 * (q2*q3 - q1*q4);
A(3,3) = q4^2 - q1^2 - q2^2 + q3^2;

end
