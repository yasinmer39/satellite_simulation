% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Appendix D - Quaternions
% Denklem D-6: Quaternion Multiplication

% Bu fonksiyon, iki quaternion'u çarpar.
% q = q1 ⊗ q2 (Hamilton çarpımı)

% q1, q2     - quaternionlar [q1, q2, q3, q4]' (q4 skaler kısım) [4x1]
% q          - çarpım quaternion [4x1]

% Quaternion çarpımı (Wertz Eq. 12-15, D-6):
%   Eğer q1 ve q2 ardışık dönüşümleri temsil ediyorsa,
%   q = q2 ⊗ q1 toplam dönüşümü verir.

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function q = quat_multiply(q1, q2)

% Vektörleri sütun vektörüne çevir
if isrow(q1)
    q1 = q1';
end
if isrow(q2)
    q2 = q2';
end

% Quaternion bileşenleri
a1 = q1(1); b1 = q1(2); c1 = q1(3); d1 = q1(4);
a2 = q2(1); b2 = q2(2); c2 = q2(3); d2 = q2(4);

% Hamilton çarpımı
q = zeros(4, 1);

q(1) = d1*a2 + a1*d2 + b1*c2 - c1*b2;
q(2) = d1*b2 - a1*c2 + b1*d2 + c1*a2;
q(3) = d1*c2 + a1*b2 - b1*a2 + c1*d2;
q(4) = d1*d2 - a1*a2 - b1*b2 - c1*c2;

% Normalize et
q = q / norm(q);

end
