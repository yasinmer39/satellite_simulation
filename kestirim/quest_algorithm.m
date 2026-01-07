% QUEST Algorithm - QUaternion ESTimator
% Wertz Chapter 12.2.3 "q Method" (Eq. 12-52 to 12-65) referans alınmıştır.
% Davenport (1968), Wahba (1965), Shuster & Oh (1981) çalışmaları temel alınmıştır.
%
% QUEST, birden fazla vektör ölçümünden optimal attitude quaternion hesaplar.
% Wahba'nın loss function'ını minimize eden quaternion'u bulur.

% v_body  - body frame'deki vektörler (normalize) [3xN]
% v_ref   - referans frame'deki vektörler (ECI) (normalize) [3xN]
% weights - her vektörün ağırlığı (opsiyonel) [1xN veya Nx1]
% q       - optimal attitude quaternion [q1,q2,q3,q4] (scalar-last) [4x1]
% A       - Direction Cosine Matrix (body <- reference) [3x3]
% loss    - Wahba loss function değeri (fit kalitesi)

% Wahba Loss Function (Wertz Eq. 12-52):
%   L(A) = Σ w_i * |v_body_i - A * v_ref_i|^2
%
% QUEST, K matrisinin en büyük özdeğerine karşılık gelen özvektörü bulur.
% K matrisi (Wertz Eq. 12-61):
%   K = [S - σI,    Z  ]
%       [Z',       σ   ]
%
% Burada:
%   B = Σ w_i * v_body_i * v_ref_i'  (Eq. 12-62a)
%   S = B + B'                        (Eq. 12-62b)
%   σ = trace(B)                      (Eq. 12-62c)
%   Z = [B23-B32, B31-B13, B12-B21]'  (Eq. 12-62d)

% Referanslar:
%   - Wertz, "Spacecraft Attitude Determination and Control", Chapter 12.2.3
%   - Shuster & Oh (1981), "Three-Axis Attitude Determination from Vector Observations"
%   - Markley & Mortari (2000), "Quaternion Attitude Estimation Using Vector Observations"

function [q, A, loss] = quest_algorithm(v_body, v_ref, weights)

% Giriş kontrolü
if nargin < 2
    error('QUEST requires at least v_body and v_ref inputs');
end

% Boyut kontrolü
[m1, n1] = size(v_body);
[m2, n2] = size(v_ref);

if m1 ~= 3 || m2 ~= 3
    error('Vectors must be 3xN matrices');
end

if n1 ~= n2
    error('Number of body and reference vectors must match');
end

num_vectors = n1;

% Varsayılan ağırlıklar (eşit)
if nargin < 3 || isempty(weights)
    weights = ones(1, num_vectors) / num_vectors;
else
    weights = weights(:)';  % Row vector
    weights = weights / sum(weights);  % Normalize
end

% Vektörleri normalize et
for i = 1:num_vectors
    v_body(:,i) = v_body(:,i) / norm(v_body(:,i));
    v_ref(:,i) = v_ref(:,i) / norm(v_ref(:,i));
end

% B matrisi hesapla (Wertz Eq. 12-62a)
B = zeros(3,3);
for i = 1:num_vectors
    B = B + weights(i) * v_body(:,i) * v_ref(:,i)';
end

% S matrisi (Wertz Eq. 12-62b)
S = B + B';

% sigma (Wertz Eq. 12-62c)
sigma = trace(B);

% Z vektörü (Wertz Eq. 12-62d)
Z = [B(2,3) - B(3,2);
     B(3,1) - B(1,3);
     B(1,2) - B(2,1)];

% K matrisi (Wertz Eq. 12-61)
K = [S - sigma*eye(3), Z;
     Z',               sigma];

% En büyük özdeğer ve özvektörü bul (Wertz Eq. 12-64, 12-65)
[V, D] = eig(K);
eigenvalues = diag(D);
[max_eigenvalue, max_idx] = max(eigenvalues);

% Optimal quaternion
q = V(:, max_idx);

% Quaternion normalizasyonu
q = q / norm(q);

% Scalar part (q4) pozitif olacak şekilde işaret ayarla
if q(4) < 0
    q = -q;
end

% Loss function değeri (Wertz Eq. 12-65)
% L = 1 - max_eigenvalue (when weights sum to 1)
loss = 1 - max_eigenvalue;

% DCM hesapla (quaternion'dan)
A = quaternion_to_dcm(q);

end


function A = quaternion_to_dcm(q)
% Quaternion'dan Direction Cosine Matrix'e dönüşüm
% Wertz Eq. 12-13a referans alınmıştır.
%
% q = [q1, q2, q3, q4] where q4 is scalar part

q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);  % scalar part

% DCM (Wertz Eq. 12-13a)
A = [q1^2 - q2^2 - q3^2 + q4^2,   2*(q1*q2 + q3*q4),           2*(q1*q3 - q2*q4);
     2*(q1*q2 - q3*q4),           -q1^2 + q2^2 - q3^2 + q4^2,  2*(q2*q3 + q1*q4);
     2*(q1*q3 + q2*q4),           2*(q2*q3 - q1*q4),           -q1^2 - q2^2 + q3^2 + q4^2];

end
