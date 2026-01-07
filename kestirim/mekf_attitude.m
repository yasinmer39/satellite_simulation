% Multiplicative Extended Kalman Filter (MEKF) for Attitude Estimation
% 
% IMU (gyro + accelerometer) ve Magnetometer füzyonu ile yönelim kestirimi.
% Quaternion hata temsili kullanılır (multiplicative error quaternion).
%
% Referanslar:
%   - Wertz Chapter 12, 17 - Attitude Determination and Prediction
%   - Markley (2003), "Attitude Error Representations for Kalman Filtering"
%   - Crassidis & Junkins (2012), "Optimal Estimation of Dynamic Systems"
%   - Lefferts et al. (1982), "Kalman Filtering for Spacecraft Attitude Estimation"

% State Vector (7x1 represented as 6x1 error state):
%   x = [δθ(3x1); b_gyro(3x1)]
%   δθ: attitude error (Gibbs vector / Rodrigues parameters)
%   b_gyro: gyro bias
%
% True State:
%   q_true = δq ⊗ q_est  (quaternion multiplication)
%   ω_true = ω_meas - b_gyro

% Fonksiyon çağrısı:
%   [q_est, state, P] = mekf_attitude(gyro, accel, mag, B_ref, g_ref, dt, state, P, params)
%
% Girdiler:
%   gyro    - gyro ölçümü [ωx, ωy, ωz] (rad/s) [3x1]
%   accel   - accelerometer ölçümü [ax, ay, az] (m/s^2) [3x1]
%   mag     - magnetometer ölçümü [Bx, By, Bz] (body frame) [3x1]
%   B_ref   - referans manyetik alan (ECI frame, normalize) [3x1]
%   g_ref   - referans gravity vektörü (ECI frame, normalize) [3x1]
%   dt      - zaman adımı (s)
%   state   - filtre durumu (struct)
%   P       - kovaryans matrisi [6x6]
%   params  - filtre parametreleri (struct)
%
% Çıktılar:
%   q_est   - kestirilen quaternion [q1,q2,q3,q4] [4x1]
%   state   - güncellenmiş filtre durumu
%   P       - güncellenmiş kovaryans matrisi

function [q_est, state, P] = mekf_attitude(gyro, accel, mag, B_ref, g_ref, dt, state, P, params)

% Varsayılan parametreler
if nargin < 9 || isempty(params)
    params = mekf_default_params();
end

% Durumu başlat (ilk çağrıda)
if nargin < 7 || isempty(state)
    state = mekf_init_state(accel, mag, B_ref, g_ref);
end

% Kovaryansı başlat (ilk çağrıda)
if nargin < 8 || isempty(P)
    P = mekf_init_covariance(params);
end

% Girdileri sütun vektörüne çevir
gyro = gyro(:);
accel = accel(:);
mag = mag(:);
B_ref = B_ref(:);
g_ref = g_ref(:);

%% ==================== PREDICTION STEP ====================
% Gyro ölçümü ile quaternion propagasyonu

% Bias-corrected angular velocity
omega = gyro - state.gyro_bias;

% Quaternion kinematic equation (Wertz Eq. 17-2)
% q_dot = 0.5 * Omega(ω) * q
q_pred = quaternion_propagate(state.q, omega, dt);

% State transition matrix F (6x6)
% Attitude error dynamics + gyro bias dynamics
F = compute_F_matrix(omega, dt, params);

% Process noise covariance Q
Q = compute_Q_matrix(dt, params);

% Covariance prediction
P_pred = F * P * F' + Q;

%% ==================== UPDATE STEP ====================
% Magnetometer ve Accelerometer ölçümleri ile güncelleme

% Normalize ölçümler
accel_norm = accel / norm(accel);
mag_norm = mag / norm(mag);
B_ref_norm = B_ref / norm(B_ref);
g_ref_norm = g_ref / norm(g_ref);

% Predicted measurements (body frame'de beklenen vektörler)
R_pred = quaternion_to_dcm(q_pred);  % ECI -> Body DCM
g_pred = R_pred * g_ref_norm;        % Predicted gravity in body
B_pred = R_pred * B_ref_norm;        % Predicted mag field in body

% Measurement residual (innovation)
% z = [accel_meas - accel_pred; mag_meas - mag_pred]
z = [accel_norm - g_pred;
     mag_norm - B_pred];

% Measurement matrix H (6x6)
H = compute_H_matrix(g_pred, B_pred);

% Measurement noise covariance R (6x6)
R = compute_R_matrix(params);

% Kalman gain
S = H * P_pred * H' + R;
K = P_pred * H' / S;

% State update
dx = K * z;  % [δθ(3x1); δb(3x1)]

% Quaternion update (multiplicative)
delta_theta = dx(1:3);
delta_q = error_to_quaternion(delta_theta);
q_est = quaternion_multiply(delta_q, q_pred);
q_est = q_est / norm(q_est);  % Normalize

% Gyro bias update
state.gyro_bias = state.gyro_bias + dx(4:6);

% Covariance update (Joseph form for numerical stability)
I_KH = eye(6) - K * H;
P = I_KH * P_pred * I_KH' + K * R * K';

% Durumu güncelle
state.q = q_est;

end


%% ==================== HELPER FUNCTIONS ====================

function state = mekf_init_state(accel, mag, B_ref, g_ref)
% TRIAD ile başlangıç yönelimi hesapla
    accel = accel(:) / norm(accel);
    mag = mag(:) / norm(mag);
    B_ref = B_ref(:) / norm(B_ref);
    g_ref = g_ref(:) / norm(g_ref);
    
    % TRIAD: gravity daha doğru, mag daha az doğru
    [q_init, ~] = triad_algorithm(accel, mag, g_ref, B_ref);
    
    state.q = q_init;
    state.gyro_bias = zeros(3,1);
end


function P = mekf_init_covariance(params)
% Başlangıç kovaryans matrisi
    P = diag([params.init_attitude_var * ones(1,3), ...
              params.init_bias_var * ones(1,3)]);
end


function params = mekf_default_params()
% MEKF varsayılan parametreleri
    
    % Process noise (gyro)
    params.gyro_noise = 0.15 * (pi/180) / sqrt(3600);  % ARW: 0.15 deg/sqrt(h)
    params.gyro_bias_noise = 0.5 * (pi/180) / 3600;    % Bias instability: 0.5 deg/h
    
    % Measurement noise
    params.accel_noise = 0.01;   % g (normalized) - çok hassas değil
    params.mag_noise = 0.05;     % normalized - daha gürültülü
    
    % Initial covariance
    params.init_attitude_var = (10 * pi/180)^2;  % 10 deg initial uncertainty
    params.init_bias_var = (1 * pi/180 / 3600)^2;  % 1 deg/h initial bias uncertainty
end


function F = compute_F_matrix(omega, dt, params)
% State transition matrix (discrete time)
% F = exp(F_c * dt) ≈ I + F_c * dt for small dt
    
    % Skew-symmetric matrix of omega
    omega_cross = skew(omega);
    
    % Continuous time F matrix
    % [δθ_dot]   [-ω×   -I] [δθ  ]
    % [b_dot ] = [0      0] [b   ]
    F_c = [-omega_cross, -eye(3);
           zeros(3,3),   zeros(3,3)];
    
    % Discrete time approximation
    F = eye(6) + F_c * dt;
end


function Q = compute_Q_matrix(dt, params)
% Process noise covariance matrix
    
    % Gyro noise affects attitude
    Q_theta = (params.gyro_noise^2 * dt) * eye(3);
    
    % Gyro bias random walk
    Q_bias = (params.gyro_bias_noise^2 * dt) * eye(3);
    
    Q = blkdiag(Q_theta, Q_bias);
end


function H = compute_H_matrix(g_pred, B_pred)
% Measurement Jacobian matrix
% h(x) = R(q) * v_ref
% dh/dδθ = -[R*v_ref]× = -[v_body_pred]×
    
    H_accel = -skew(g_pred);
    H_mag = -skew(B_pred);
    
    % Bias'ın ölçüm üzerinde direkt etkisi yok
    H = [H_accel, zeros(3,3);
         H_mag,   zeros(3,3)];
end


function R = compute_R_matrix(params)
% Measurement noise covariance matrix
    R_accel = params.accel_noise^2 * eye(3);
    R_mag = params.mag_noise^2 * eye(3);
    
    R = blkdiag(R_accel, R_mag);
end


function q_new = quaternion_propagate(q, omega, dt)
% Quaternion propagation using angular velocity
% q_dot = 0.5 * Omega(ω) * q
    
    % Omega matrix
    wx = omega(1); wy = omega(2); wz = omega(3);
    Omega = [0,   wz, -wy,  wx;
            -wz,  0,   wx,  wy;
             wy, -wx,  0,   wz;
            -wx, -wy, -wz,  0];
    
    % First-order integration (good for small dt)
    q_dot = 0.5 * Omega * q;
    q_new = q + q_dot * dt;
    q_new = q_new / norm(q_new);
end


function A = quaternion_to_dcm(q)
% Quaternion'dan DCM'e dönüşüm
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    A = [q1^2-q2^2-q3^2+q4^2,   2*(q1*q2+q3*q4),       2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4),       -q1^2+q2^2-q3^2+q4^2,  2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4),       2*(q2*q3-q1*q4),       -q1^2-q2^2+q3^2+q4^2];
end


function q = error_to_quaternion(delta_theta)
% Küçük açı hatasından quaternion'a dönüşüm
% δq ≈ [δθ/2; 1] for small δθ
    
    delta_theta = delta_theta(:);
    norm_dt = norm(delta_theta);
    
    if norm_dt < 1e-10
        q = [0; 0; 0; 1];
    else
        q = [delta_theta/2; 1];
        q = q / norm(q);
    end
end


function q = quaternion_multiply(p, q)
% Quaternion multiplication: r = p ⊗ q
    p1 = p(1); p2 = p(2); p3 = p(3); p4 = p(4);
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    r = [p4*q1 + p1*q4 + p2*q3 - p3*q2;
         p4*q2 - p1*q3 + p2*q4 + p3*q1;
         p4*q3 + p1*q2 - p2*q1 + p3*q4;
         p4*q4 - p1*q1 - p2*q2 - p3*q3];
    q = r;
end


function S = skew(v)
% Skew-symmetric (cross product) matrix
    S = [0,    -v(3),  v(2);
         v(3),  0,    -v(1);
        -v(2),  v(1),  0];
end
