% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Birleşik Yönelim Dinamiği ODE Fonksiyonu

% Bu fonksiyon, yönelim dinamiği için state-space ODE formunu sağlar.
% MATLAB ode45/ode113 gibi integratörlerle kullanılabilir.

% t         - zaman (s) [skaler]
% state     - durum vektörü [7x1]:
%             [q1, q2, q3, q4, omega1, omega2, omega3]'
%             q: quaternion (birim)
%             omega: açısal hız body frame'de (rad/s)
% I         - asal atalet momentleri [I1, I2, I3] (kg·m²) [3x1]
% N_func    - tork fonksiyonu handle: N = N_func(t, state)
%             veya sabit tork vektörü [3x1]
% state_dot - durum türevi [7x1]

% Denklemler:
%   Kinematik: dq/dt = (1/2) * Ω * q        (Wertz Eq. 16-3)
%   Dinamik:   I*dω/dt = N - ω × (I*ω)      (Wertz Eq. 16-48)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar:
%   quaternion_kinematics, euler_equations

function state_dot = attitude_dynamics_ode(t, state, I, N_func)

% State vektörünü ayır
q = state(1:4);
omega = state(5:7);

% Quaternion normalizasyonu (sayısal drift'i önlemek için)
q = q / norm(q);

% Tork hesaplama
if isa(N_func, 'function_handle')
    N = N_func(t, state);
else
    N = N_func;  % Sabit tork
end

% Kinematik denklemler (quaternion)
q_dot = quaternion_kinematics(q, omega);

% Dinamik denklemler (Euler)
omega_dot = euler_equations(omega, I, N);

% State türevi
state_dot = [q_dot; omega_dot];

end
