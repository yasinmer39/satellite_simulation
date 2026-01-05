% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Appendix D - Quaternions
% Denklem D-3: Quaternion Conjugate/Inverse

% Bu fonksiyon, quaternion'un tersini (conjugate) hesaplar.
% Birim quaternion için: q* = q^(-1)

% q         - quaternion [q1, q2, q3, q4]' [4x1]
% q_inv     - ters quaternion [4x1]

% Quaternion conjugate (Wertz Eq. D-3):
%   q* = [q4 - i*q1 - j*q2 - k*q3]
%      = [-q1, -q2, -q3, q4]'

% Özellik: q ⊗ q* = [0, 0, 0, 1]' (identity)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function q_inv = quat_inverse(q)

% Vektörü sütun vektörüne çevir
if isrow(q)
    q = q';
end

% Conjugate (birim quaternion için ters)
q_inv = [-q(1); -q(2); -q(3); q(4)];

% Birim quaternion değilse normalize et
q_inv = q_inv / (norm(q)^2);

end
