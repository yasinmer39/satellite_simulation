% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Denklem 16-7: Direction Cosine Matrix Kinematic Equations

% Bu fonksiyon, DCM (Direction Cosine Matrix) kinematik denklemlerini hesaplar.

% A         - DCM (yön kosinüsleri matrisi) [3x3]
%             Body frame'den Reference frame'e dönüşüm
% omega     - açısal hız vektörü body frame'de (rad/s) [3x1]
% A_dot     - DCM türevi [3x3]

% DCM Kinematik (Wertz Eq. 16-7):
%   dA/dt = Ω' * A
%
% Burada Ω', 3x3 skew-symmetric matristir (Wertz Eq. 16-6b):
%       [  0   -ω3   ω2  ]
%   Ω' = [ ω3    0   -ω1  ]
%       [ -ω2   ω1    0   ]

% Not: Bu denklem body frame'deki açısal hızı kullanır.
%      Reference frame'deki açısal hız için: dA/dt = A * Ω

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: skew_symmetric

function A_dot = dcm_kinematics(A, omega)

% Vektörü sütun vektörüne çevir
if isrow(omega)
    omega = omega';
end

% Omega' skew-symmetric matrisi (Wertz Eq. 16-6b)
Omega_prime = skew_symmetric(omega);

% DCM türevi (Wertz Eq. 16-7)
A_dot = Omega_prime * A;

end
