% Spacecraft Attitude Determination and Control (Wertz) kitabından uyarlanmıştır.
% Appendix D - Quaternions

% Bu fonksiyon, quaternion'u normalize eder.
% |q| = 1 olacak şekilde.

% q         - quaternion [q1, q2, q3, q4]' [4x1]
% q_norm    - normalize edilmiş quaternion [4x1]

% Birim quaternion özelliği:
%   |q|² = q1² + q2² + q3² + q4² = 1

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function q_norm = quat_normalize(q)

% Vektörü sütun vektörüne çevir
if isrow(q)
    q = q';
end

% Normalize et
q_mag = norm(q);

if q_mag < 1e-10
    warning('Quaternion magnitude çok küçük!');
    q_norm = [0; 0; 0; 1];  % identity quaternion
else
    q_norm = q / q_mag;
end

% q4 pozitif olsun (opsiyonel convention)
% if q_norm(4) < 0
%     q_norm = -q_norm;
% end

end
