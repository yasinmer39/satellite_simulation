% TRIAD Algorithm - Three-Axis Attitude Determination
% Wertz Chapter 12.2.2 Algebraic Method (Eq. 12-39 to 12-45) referans alınmıştır.
% 
% İki vektör ölçümünden attitude matrisi (DCM) ve quaternion hesaplar.
% TRIAD, iki referans vektör ve bunların body frame'deki ölçümlerini
% kullanarak spacecraft yönelimini belirler.

% v1_body  - birinci vektör body frame'de (normalize) [3x1]
% v2_body  - ikinci vektör body frame'de (normalize) [3x1]
% v1_ref   - birinci vektör referans frame'de (ECI) (normalize) [3x1]
% v2_ref   - ikinci vektör referans frame'de (ECI) (normalize) [3x1]
% q        - attitude quaternion [q1,q2,q3,q4] (scalar-last) [4x1]
% A        - Direction Cosine Matrix (body <- reference) [3x3]

% TRIAD Algoritması (Wertz Eq. 12-39 to 12-45):
%   1. İki vektörden ortonormal baz oluştur (body ve reference)
%   2. M_body = [t1_b, t2_b, t3_b]
%   3. M_ref  = [t1_r, t2_r, t3_r]
%   4. A = M_body * M_ref'
%
% Birinci vektör (v1) daha doğru kabul edilir, ikinci vektör (v2)
% sadece v1 etrafındaki rotasyonu belirler.

% Not: Magnetometer daha az doğru olduğundan v2 olarak kullanılmalı.
%      IMU accelerometer gravity vektörü v1 olarak kullanılabilir.

% Referanslar:
%   - Wertz, "Spacecraft Attitude Determination and Control", Chapter 12.2.2
%   - Shuster & Oh (1981), "Three-Axis Attitude Determination from Vector Observations"

function [q, A] = triad_algorithm(v1_body, v2_body, v1_ref, v2_ref)

% Giriş kontrolü
if nargin < 4
    error('TRIAD requires 4 input vectors: v1_body, v2_body, v1_ref, v2_ref');
end

% Sütun vektörlerine çevir
v1_body = v1_body(:);
v2_body = v2_body(:);
v1_ref = v1_ref(:);
v2_ref = v2_ref(:);

% Vektörleri normalize et
v1_body = v1_body / norm(v1_body);
v2_body = v2_body / norm(v2_body);
v1_ref = v1_ref / norm(v1_ref);
v2_ref = v2_ref / norm(v2_ref);

% Vektörlerin paralel olup olmadığını kontrol et (Wertz Eq. 12-40)
dot_body = abs(dot(v1_body, v2_body));
dot_ref = abs(dot(v1_ref, v2_ref));

if dot_body > 0.999 || dot_ref > 0.999
    warning('TRIAD: Vectors are nearly parallel, attitude determination may be inaccurate');
end

% TRIAD ortonormal baz oluşturma (Wertz Eq. 12-39a,b,c)
% Body frame'de baz vektörleri
t1_body = v1_body;                                    % Eq. 12-39a
t2_body = cross(v1_body, v2_body);                    % Eq. 12-39b (normalize edilmemiş)
t2_body = t2_body / norm(t2_body);                    % Normalize
t3_body = cross(t1_body, t2_body);                    % Eq. 12-39c

% Reference frame'de baz vektörleri
t1_ref = v1_ref;
t2_ref = cross(v1_ref, v2_ref);
t2_ref = t2_ref / norm(t2_ref);
t3_ref = cross(t1_ref, t2_ref);

% Body ve Reference matrisleri (Wertz Eq. 12-41, 12-42)
M_body = [t1_body, t2_body, t3_body];
M_ref = [t1_ref, t2_ref, t3_ref];

% Attitude matrix (DCM) hesaplama (Wertz Eq. 12-45)
% A: body <- reference dönüşümü
% v_body = A * v_ref
A = M_body * M_ref';

% DCM'in ortogonalliğini kontrol et ve düzelt
% SVD ile en yakın ortogonal matrise projeksiyon
[U, ~, V] = svd(A);
A = U * V';

% Determinant kontrolü (proper rotation için +1 olmalı)
if det(A) < 0
    A = -A;
    warning('TRIAD: DCM determinant was negative, sign corrected');
end

% DCM'den quaternion'a dönüşüm (Shuster method)
q = dcm_to_quaternion(A);

end


function q = dcm_to_quaternion(A)
% Direction Cosine Matrix'ten Quaternion'a dönüşüm
% Wertz Eq. 12-14a,b,c,d ve Shuster (1993) referans alınmıştır.
% 
% Quaternion convention: q = [q1, q2, q3, q4] where q4 is scalar part
% q4 = cos(theta/2), [q1,q2,q3] = e*sin(theta/2)

% Trace hesapla
tr = trace(A);

% En büyük quaternion bileşenini bul (sayısal kararlılık için)
% Wertz'deki 4 alternatif formülden en kararlı olanı seç
K = [A(1,1)-A(2,2)-A(3,3),  A(2,1)+A(1,2),          A(3,1)+A(1,3),          A(2,3)-A(3,2);
     A(2,1)+A(1,2),          A(2,2)-A(1,1)-A(3,3),  A(3,2)+A(2,3),          A(3,1)-A(1,3);
     A(3,1)+A(1,3),          A(3,2)+A(2,3),          A(3,3)-A(1,1)-A(2,2),  A(1,2)-A(2,1);
     A(2,3)-A(3,2),          A(3,1)-A(1,3),          A(1,2)-A(2,1),          tr];

K = K / 3;

% En büyük özdeğere karşılık gelen özvektör = quaternion
[V, D] = eig(K);
[~, max_idx] = max(diag(D));
q = V(:, max_idx);

% Quaternion normalizasyonu
q = q / norm(q);

% Scalar part (q4) pozitif olacak şekilde işaret ayarla
if q(4) < 0
    q = -q;
end

end
