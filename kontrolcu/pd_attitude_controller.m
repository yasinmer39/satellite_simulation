% PD Attitude Controller - Quaternion Feedback
% Case Study (g): PID kontrolcü tasarlayınız
%
% Bu fonksiyon, quaternion hata geri beslemesi kullanarak
% kontrol torku hesaplar.
%
% Referanslar:
%   - Wertz Chapter 18, Eq. 18-8: Nc = -K(τθ̇ + θ)
%   - Wie & Barba (1985): "Quaternion Feedback for Spacecraft Large Angle Maneuvers"
%   - Sidi, "Spacecraft Dynamics and Control"
%
% Kontrol Kanunu:
%   τ_cmd = -Kp * q_err_v - Kd * ω_err + ω × (I*ω + h_w)
%
% Burada:
%   q_err_v = quaternion error'un vektör kısmı [3x1]
%   ω_err   = angular velocity error [3x1]
%   Kp, Kd  = proportional ve derivative gain matrisleri [3x3]
%
% Gain Tasarımı (Wertz Eq. 18-12):
%   ωn = sqrt(Kp/I)   → Natural frequency
%   ζ  = Kd/(2*sqrt(Kp*I)) → Damping ratio
%
% Author: Mehmet Yasin Meriç
% Date: 2026
% Course: Uzay Sistemleri Vaka Çalışması

function [tau_cmd, ctrl_info] = pd_attitude_controller(q_current, q_desired, omega_current, omega_desired, params)
% PD_ATTITUDE_CONTROLLER - Quaternion-based PD kontrol
%
% Girdiler:
%   q_current    - Mevcut quaternion (ECI→Body) [4x1]
%   q_desired    - Hedef quaternion (ECI→Body) [4x1]
%   omega_current- Mevcut açısal hız (body frame, rad/s) [3x1]
%   omega_desired- Hedef açısal hız (body frame, rad/s) [3x1]
%   params       - Kontrolcü parametreleri (struct):
%                  .Kp     - Proportional gain [3x3] veya [scalar]
%                  .Kd     - Derivative gain [3x3] veya [scalar]
%                  .I      - Inertia matrix [3x3]
%                  .h_wheels - Wheel angular momentum [3x1] (optional)
%                  .tau_max - Maximum torque [3x1] veya [scalar]
%                  .feedforward - Enable feedforward compensation (bool)
%
% Çıktılar:
%   tau_cmd   - Komut edilen kontrol torku (body frame, Nm) [3x1]
%   ctrl_info - Kontrol bilgileri (debug için) struct
%
% Quaternion convention: q = [q1; q2; q3; q4] (vector-first, scalar-last)

    %% Input validation and defaults
    if nargin < 5
        error('All inputs required: q_current, q_desired, omega_current, omega_desired, params');
    end
    
    % Ensure column vectors
    q_current = q_current(:);
    q_desired = q_desired(:);
    omega_current = omega_current(:);
    omega_desired = omega_desired(:);
    
    % Normalize quaternions
    q_current = q_current / norm(q_current);
    q_desired = q_desired / norm(q_desired);
    
    % Extract parameters
    I = params.I;
    
    % Handle scalar or matrix gains
    if isscalar(params.Kp)
        Kp = params.Kp * eye(3);
    else
        Kp = params.Kp;
    end
    
    if isscalar(params.Kd)
        Kd = params.Kd * eye(3);
    else
        Kd = params.Kd;
    end
    
    % Optional parameters
    if isfield(params, 'h_wheels')
        h_wheels = params.h_wheels(:);
    else
        h_wheels = [0; 0; 0];
    end
    
    if isfield(params, 'tau_max')
        if isscalar(params.tau_max)
            tau_max = params.tau_max * ones(3, 1);
        else
            tau_max = params.tau_max(:);
        end
    else
        tau_max = inf * ones(3, 1);
    end
    
    if isfield(params, 'feedforward')
        use_feedforward = params.feedforward;
    else
        use_feedforward = true;
    end
    
    %% Compute quaternion error
    % q_error = q_desired^(-1) ⊗ q_current
    % Bu, body frame'de hatayı verir
    %
    % Wertz ve Wie'ye göre:
    % q_err = q_d* ⊗ q (multiplicative error)
    
    q_desired_conj = [-q_desired(1:3); q_desired(4)];  % Quaternion conjugate
    q_error = quat_multiply(q_desired_conj, q_current);
    
    % Ensure shortest path (avoid unwinding)
    % If scalar part is negative, flip sign
    if q_error(4) < 0
        q_error = -q_error;
    end
    
    % Extract vector part of quaternion error
    q_error_v = q_error(1:3);  % This is ≈ 0.5 * theta_error for small angles
    q_error_s = q_error(4);    % Scalar part (≈ 1 for small angles)
    
    % Attitude error in radians (small angle approximation)
    % theta_error ≈ 2 * q_error_v (for small errors)
    % For large angles: theta = 2 * acos(q_error_s)
    theta_error_mag = 2 * acos(min(1, abs(q_error_s)));
    
    %% Compute angular velocity error
    omega_error = omega_current - omega_desired;
    
    %% PD Control Law
    % Wie (1985), Sidi (1997):
    % τ = -Kp * sign(q4) * q_v - Kd * ω_error
    %
    % Full form with feedforward (computed torque):
    % τ = -Kp * q_v - Kd * ω_error + ω × (I*ω + h_w)
    
    % Proportional term (attitude error)
    % Note: sign(q_error_s) ensures shortest rotation path
    tau_p = -Kp * sign(q_error_s) * q_error_v;
    
    % Derivative term (rate error)
    tau_d = -Kd * omega_error;
    
    % Feedforward term (gyroscopic compensation)
    if use_feedforward
        L_total = I * omega_current + h_wheels;  % Total angular momentum
        tau_ff = cross(omega_current, L_total);
    else
        tau_ff = [0; 0; 0];
    end
    
    % Total control torque
    tau_cmd = tau_p + tau_d + tau_ff;
    
    %% Apply saturation
    for i = 1:3
        if abs(tau_cmd(i)) > tau_max(i)
            tau_cmd(i) = sign(tau_cmd(i)) * tau_max(i);
        end
    end
    
    %% Output info
    ctrl_info.q_error = q_error;
    ctrl_info.q_error_v = q_error_v;
    ctrl_info.theta_error_deg = theta_error_mag * 180/pi;
    ctrl_info.omega_error = omega_error;
    ctrl_info.tau_p = tau_p;
    ctrl_info.tau_d = tau_d;
    ctrl_info.tau_ff = tau_ff;
    ctrl_info.saturated = any(abs(tau_cmd) >= tau_max - 1e-10);
    
end

%% ==================== QUATERNION MULTIPLICATION ====================
function r = quat_multiply(p, q)
% Quaternion multiplication: r = p ⊗ q
% Convention: q = [q1; q2; q3; q4] = [qv; qs] (scalar-last)
%
% Reference: Wertz Appendix E

    p1 = p(1); p2 = p(2); p3 = p(3); p4 = p(4);
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    r = [p4*q1 + p1*q4 + p2*q3 - p3*q2;
         p4*q2 - p1*q3 + p2*q4 + p3*q1;
         p4*q3 + p1*q2 - p2*q1 + p3*q4;
         p4*q4 - p1*q1 - p2*q2 - p3*q3];
end

%% ==================== GAIN CALCULATION HELPER ====================
function [Kp, Kd] = compute_pd_gains(I, zeta, omega_n)
% COMPUTE_PD_GAINS - Hesapla PD gains from damping ratio ve natural frequency
%
% Wertz Eq. 18-12:
%   ωn = sqrt(K/I)
%   ζ  = Kτ/(2*sqrt(K*I))
%
% Çözüm:
%   Kp = I * ωn²
%   Kd = 2 * ζ * ωn * I
%
% Girdiler:
%   I       - Moment of inertia (scalar veya [3x3])
%   zeta    - Damping ratio (tipik: 0.7-1.0)
%   omega_n - Natural frequency (rad/s)
%
% Çıktılar:
%   Kp - Proportional gain
%   Kd - Derivative gain

    if isscalar(I)
        Kp = I * omega_n^2;
        Kd = 2 * zeta * omega_n * I;
    else
        % For matrix I, use diagonal terms
        Kp = diag(diag(I)) * omega_n^2;
        Kd = 2 * zeta * omega_n * diag(diag(I));
    end
end
