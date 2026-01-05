% Spacecraft Attitude Determination and Control (Wertz) kitabından uyarlanmıştır.
% Section 12.1 - Parameterization of the Attitude

% Bu fonksiyon, Direction Cosine Matrix (DCM)'den quaternion hesaplar.
% Shepperd's method kullanılır (sayısal kararlılık için).

% A         - DCM (yön kosinüsleri matrisi) [3x3]
% q         - quaternion [q1, q2, q3, q4]' (q4 skaler kısım) [4x1]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function q = dcm2quat_wertz(A)

% Trace hesapla
tr = trace(A);

% Shepperd's method: en büyük elemanı bul
% Bu sayısal kararlılık sağlar

% Dört olası değer
K = zeros(4, 1);
K(1) = 1 + 2*A(1,1) - tr;
K(2) = 1 + 2*A(2,2) - tr;
K(3) = 1 + 2*A(3,3) - tr;
K(4) = 1 + tr;

% En büyük değerin indeksini bul
[~, max_idx] = max(K);

q = zeros(4, 1);

switch max_idx
    case 1
        q(1) = sqrt(K(1)) / 2;
        q(2) = (A(1,2) + A(2,1)) / (4 * q(1));
        q(3) = (A(1,3) + A(3,1)) / (4 * q(1));
        q(4) = (A(2,3) - A(3,2)) / (4 * q(1));
        
    case 2
        q(2) = sqrt(K(2)) / 2;
        q(1) = (A(1,2) + A(2,1)) / (4 * q(2));
        q(3) = (A(2,3) + A(3,2)) / (4 * q(2));
        q(4) = (A(3,1) - A(1,3)) / (4 * q(2));
        
    case 3
        q(3) = sqrt(K(3)) / 2;
        q(1) = (A(1,3) + A(3,1)) / (4 * q(3));
        q(2) = (A(2,3) + A(3,2)) / (4 * q(3));
        q(4) = (A(1,2) - A(2,1)) / (4 * q(3));
        
    case 4
        q(4) = sqrt(K(4)) / 2;
        q(1) = (A(2,3) - A(3,2)) / (4 * q(4));
        q(2) = (A(3,1) - A(1,3)) / (4 * q(4));
        q(3) = (A(1,2) - A(2,1)) / (4 * q(4));
end

% Normalize et
q = q / norm(q);

% q4 pozitif olsun (convention)
if q(4) < 0
    q = -q;
end

end
