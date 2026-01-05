% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Denklem 16-3: Quaternion Kinematic Equations

% Bu fonksiyon, quaternion kinematik denklemlerini hesaplar.
% Singularity (gimbal lock) problemi yoktur.

% q         - quaternion [q1, q2, q3, q4]' (q4 skaler kısım) [4x1]
% omega     - açısal hız vektörü body frame'de (rad/s) [3x1]
% q_dot     - quaternion türevi [4x1]

% Quaternion Kinematik (Wertz Eq. 16-3):
%   dq/dt = (1/2) * Ω * q
%
% Burada Ω, 4x4 skew-symmetric matristir (Wertz Eq. 16-2b):
%       [  0    ω3   -ω2   ω1  ]
%   Ω = [ -ω3   0     ω1   ω2  ]
%       [  ω2  -ω1    0    ω3  ]
%       [ -ω1  -ω2   -ω3   0   ]

% Quaternion Convention:
%   q = [q1, q2, q3, q4]' = [q_vec; q_scalar]
%   |q| = 1 (birim quaternion)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function q_dot = quaternion_kinematics(q, omega)

% Vektörleri sütun vektörüne çevir
if isrow(q)
    q = q';
end
if isrow(omega)
    omega = omega';
end

% Açısal hız bileşenleri
wx = omega(1);
wy = omega(2);
wz = omega(3);

% Omega matrisi (Wertz Eq. 16-2b)
Omega = [  0    wz   -wy   wx;
          -wz   0     wx   wy;
           wy  -wx    0    wz;
          -wx  -wy   -wz   0  ];

% Quaternion türevi (Wertz Eq. 16-3)
q_dot = 0.5 * Omega * q;

end
