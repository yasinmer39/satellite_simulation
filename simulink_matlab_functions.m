%% ========================================================================
%  SIMULINK MATLAB FUNCTION BLOKLARI
%  ========================================================================
%  Bu dosya, Simulink'te MATLAB Function bloklarına kopyalanacak
%  tüm fonksiyonları içerir.
%
%  Her bölümü ilgili MATLAB Function bloğuna kopyalayın.
%
%  Author: Mehmet Yasin Meriç
%  Date: 2026
%% ========================================================================

%% ========================================================================
%  1. ORBIT PROPAGATOR
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Orbit_Propagator"
%  
%  Girdiler: t (scalar), orbit (bus/struct)
%  Çıktılar: r_eci (3x1), v_eci (3x1)

function [r_eci, v_eci] = fcn(t, a, e, i, RAAN, omega_arg, M0, n)
%#codegen
% ORBIT_PROPAGATOR - Circular orbit propagation
%
% Inputs:
%   t         - Simulation time (s)
%   a         - Semi-major axis (m)
%   e         - Eccentricity (should be 0 for circular)
%   i         - Inclination (rad)
%   RAAN      - Right ascension of ascending node (rad)
%   omega_arg - Argument of perigee (rad)
%   M0        - Initial mean anomaly (rad)
%   n         - Mean motion (rad/s)

    % Mean anomaly at current time
    M = n * t + M0;
    
    % For circular orbit, true anomaly = mean anomaly
    nu = M;
    
    % Position in perifocal (orbital plane) coordinates
    r_peri = a * [cos(nu); sin(nu); 0];
    v_peri = sqrt(398600.4418e9 / a) * [-sin(nu); cos(nu); 0];
    
    % Rotation matrices (perifocal to ECI)
    % R = R3(-RAAN) * R1(-i) * R3(-omega)
    
    c_O = cos(RAAN); s_O = sin(RAAN);
    c_i = cos(i);    s_i = sin(i);
    c_w = cos(omega_arg); s_w = sin(omega_arg);
    
    % Combined rotation matrix (3-1-3 sequence)
    R11 = c_O*c_w - s_O*c_i*s_w;
    R12 = -c_O*s_w - s_O*c_i*c_w;
    R13 = s_O*s_i;
    R21 = s_O*c_w + c_O*c_i*s_w;
    R22 = -s_O*s_w + c_O*c_i*c_w;
    R23 = -c_O*s_i;
    R31 = s_i*s_w;
    R32 = s_i*c_w;
    R33 = c_i;
    
    % Transform to ECI
    r_eci = [R11*r_peri(1) + R12*r_peri(2) + R13*r_peri(3);
             R21*r_peri(1) + R22*r_peri(2) + R23*r_peri(3);
             R31*r_peri(1) + R32*r_peri(2) + R33*r_peri(3)];
    
    v_eci = [R11*v_peri(1) + R12*v_peri(2) + R13*v_peri(3);
             R21*v_peri(1) + R22*v_peri(2) + R23*v_peri(3);
             R31*v_peri(1) + R32*v_peri(2) + R33*v_peri(3)];
end


%% ========================================================================
%  2. SUN POSITION
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Sun_Position"

function sun_eci = fcn(t, jd0)
%#codegen
% SUN_POSITION - Compute sun direction in ECI frame
%
% Inputs:
%   t   - Simulation time (s)
%   jd0 - Julian date at t=0

    % Julian date at current time
    JD = jd0 + t / 86400;
    
    % Julian centuries from J2000
    T = (JD - 2451545.0) / 36525.0;
    
    % Mean longitude and anomaly (degrees)
    L0 = mod(280.46 + 36000.77 * T, 360.0);
    M = mod(357.53 + 35999.05 * T, 360.0) * 0.017453292519943;  % to rad
    
    % Ecliptic longitude (radians)
    lambda = (L0 + 1.915 * sin(M) + 0.02 * sin(2*M)) * 0.017453292519943;
    
    % Obliquity of ecliptic (radians)
    eps = (23.439 - 0.013 * T) * 0.017453292519943;
    
    % Sun direction in ECI (unit vector)
    sun_eci = [cos(lambda);
               cos(eps) * sin(lambda);
               sin(eps) * sin(lambda)];
    
    % Normalize
    sun_mag = sqrt(sun_eci(1)^2 + sun_eci(2)^2 + sun_eci(3)^2);
    sun_eci = sun_eci / sun_mag;
end


%% ========================================================================
%  3. MAGNETIC FIELD (IGRF Dipole)
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Magnetic_Field"

function B_eci = fcn(r_eci, B0, Re)
%#codegen
% MAGNETIC_FIELD - Earth magnetic field (dipole model)
%
% Inputs:
%   r_eci - Position in ECI (m)
%   B0    - Magnetic field constant (T), typically 3.12e-5
%   Re    - Earth radius (m)

    % Distance from Earth center
    r_mag = sqrt(r_eci(1)^2 + r_eci(2)^2 + r_eci(3)^2);
    
    % Unit position vector
    r_hat = r_eci / r_mag;
    
    % Magnetic dipole aligned with Z axis
    m_hat = [0; 0; 1];
    
    % Dipole field equation
    % B = B0 * (Re/r)^3 * [3*(m·r_hat)*r_hat - m]
    m_dot_r = m_hat(3) * r_hat(3);  % Since m_hat = [0;0;1]
    factor = B0 * (Re / r_mag)^3;
    
    B_eci = factor * (3 * m_dot_r * r_hat - m_hat);
end


%% ========================================================================
%  4. SPACECRAFT DYNAMICS
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Spacecraft_Dynamics"
%
%  NOT: Bu fonksiyon, Integrator bloklarından ÖNCE kullanılır.
%  Integrator'ler omega ve q'yu integrate eder.

function [omega_dot, q_dot] = fcn(omega, q, tau_total, h_wheels, I, I_inv)
%#codegen
% SPACECRAFT_DYNAMICS - Euler's equation and quaternion kinematics
%
% Inputs:
%   omega     - Angular velocity (rad/s) [3x1]
%   q         - Quaternion [4x1] scalar-last
%   tau_total - Total torque (Nm) [3x1]
%   h_wheels  - Wheel angular momentum (Nms) [3x1]
%   I         - Inertia matrix [3x3]
%   I_inv     - Inverse inertia matrix [3x3]

    % Total angular momentum
    L = I * omega + h_wheels;
    
    % Euler's equation: I*omega_dot = tau - omega x L
    omega_cross_L = [omega(2)*L(3) - omega(3)*L(2);
                     omega(3)*L(1) - omega(1)*L(3);
                     omega(1)*L(2) - omega(2)*L(1)];
    
    omega_dot = I_inv * (tau_total - omega_cross_L);
    
    % Quaternion kinematics: q_dot = 0.5 * Omega(omega) * q
    % Omega matrix
    w1 = omega(1); w2 = omega(2); w3 = omega(3);
    
    q_dot = 0.5 * [ w3*q(2) - w2*q(3) + w1*q(4);
                   -w3*q(1) + w1*q(3) + w2*q(4);
                    w2*q(1) - w1*q(2) + w3*q(4);
                   -w1*q(1) - w2*q(2) - w3*q(3)];
end


%% ========================================================================
%  5. QUATERNION NORMALIZATION
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Quat_Normalize"
%  
%  Bu bloğu Quaternion Integrator'ün çıkışına bağlayın.

function q_norm = fcn(q)
%#codegen
% QUAT_NORMALIZE - Normalize quaternion
    
    mag = sqrt(q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2);
    q_norm = q / mag;
    
    % Ensure scalar part positive (avoid unwinding)
    if q_norm(4) < 0
        q_norm = -q_norm;
    end
end


%% ========================================================================
%  6. DCM FROM QUATERNION
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Quat_to_DCM"

function A = fcn(q)
%#codegen
% QUAT_TO_DCM - Convert quaternion to Direction Cosine Matrix
%
% Input:  q = [q1; q2; q3; q4] scalar-last
% Output: A = DCM (3x3), transforms from ECI to Body

    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    A = [q4^2+q1^2-q2^2-q3^2,  2*(q1*q2+q3*q4),      2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4),      q4^2-q1^2+q2^2-q3^2,  2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4),      2*(q2*q3-q1*q4),      q4^2-q1^2-q2^2+q3^2];
end


%% ========================================================================
%  7. SUN POINTING GUIDANCE
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Sun_Pointing_Guidance"

function q_cmd = fcn(sun_eci, body_axis)
%#codegen
% SUN_POINTING_GUIDANCE - Compute desired quaternion for sun pointing
%
% Inputs:
%   sun_eci   - Sun direction in ECI [3x1]
%   body_axis - Body axis to point at sun [3x1], typically [0;0;1]

    % Normalize inputs
    sun_mag = sqrt(sun_eci(1)^2 + sun_eci(2)^2 + sun_eci(3)^2);
    sun_n = sun_eci / sun_mag;
    
    body_mag = sqrt(body_axis(1)^2 + body_axis(2)^2 + body_axis(3)^2);
    body_n = body_axis / body_mag;
    
    % Build reference frame in ECI (target)
    z_eci = sun_n;
    
    % Choose secondary axis perpendicular to sun
    if abs(z_eci(3)) < 0.99
        ref_sec = [0; 0; 1];
    else
        ref_sec = [0; 1; 0];
    end
    
    % Gram-Schmidt
    y_eci = [z_eci(2)*ref_sec(3) - z_eci(3)*ref_sec(2);
             z_eci(3)*ref_sec(1) - z_eci(1)*ref_sec(3);
             z_eci(1)*ref_sec(2) - z_eci(2)*ref_sec(1)];
    y_mag = sqrt(y_eci(1)^2 + y_eci(2)^2 + y_eci(3)^2);
    y_eci = y_eci / y_mag;
    
    x_eci = [y_eci(2)*z_eci(3) - y_eci(3)*z_eci(2);
             y_eci(3)*z_eci(1) - y_eci(1)*z_eci(3);
             y_eci(1)*z_eci(2) - y_eci(2)*z_eci(1)];
    
    % Build body frame
    z_body = body_n;
    if abs(z_body(3)) < 0.99
        ref_b = [0; 0; 1];
    else
        ref_b = [0; 1; 0];
    end
    
    y_body = [z_body(2)*ref_b(3) - z_body(3)*ref_b(2);
              z_body(3)*ref_b(1) - z_body(1)*ref_b(3);
              z_body(1)*ref_b(2) - z_body(2)*ref_b(1)];
    y_mag_b = sqrt(y_body(1)^2 + y_body(2)^2 + y_body(3)^2);
    y_body = y_body / y_mag_b;
    
    x_body = [y_body(2)*z_body(3) - y_body(3)*z_body(2);
              y_body(3)*z_body(1) - y_body(1)*z_body(3);
              y_body(1)*z_body(2) - y_body(2)*z_body(1)];
    
    % DCM: A_eci = [x_eci, y_eci, z_eci]' (ECI basis as rows)
    % DCM: A_body = [x_body, y_body, z_body] (Body basis as columns)
    % A_eci2body = A_body' * A_eci
    
    % Build DCM
    A11 = x_body(1)*x_eci(1) + y_body(1)*y_eci(1) + z_body(1)*z_eci(1);
    A12 = x_body(1)*x_eci(2) + y_body(1)*y_eci(2) + z_body(1)*z_eci(2);
    A13 = x_body(1)*x_eci(3) + y_body(1)*y_eci(3) + z_body(1)*z_eci(3);
    A21 = x_body(2)*x_eci(1) + y_body(2)*y_eci(1) + z_body(2)*z_eci(1);
    A22 = x_body(2)*x_eci(2) + y_body(2)*y_eci(2) + z_body(2)*z_eci(2);
    A23 = x_body(2)*x_eci(3) + y_body(2)*y_eci(3) + z_body(2)*z_eci(3);
    A31 = x_body(3)*x_eci(1) + y_body(3)*y_eci(1) + z_body(3)*z_eci(1);
    A32 = x_body(3)*x_eci(2) + y_body(3)*y_eci(2) + z_body(3)*z_eci(2);
    A33 = x_body(3)*x_eci(3) + y_body(3)*y_eci(3) + z_body(3)*z_eci(3);
    
    % DCM to Quaternion (Shepperd's method)
    tr = A11 + A22 + A33;
    
    S1 = 1 + tr;
    S2 = 1 + 2*A11 - tr;
    S3 = 1 + 2*A22 - tr;
    S4 = 1 + 2*A33 - tr;
    
    if S1 >= S2 && S1 >= S3 && S1 >= S4
        q4 = 0.5 * sqrt(S1);
        q1 = (A23 - A32) / (4*q4);
        q2 = (A31 - A13) / (4*q4);
        q3 = (A12 - A21) / (4*q4);
    elseif S2 >= S1 && S2 >= S3 && S2 >= S4
        q1 = 0.5 * sqrt(S2);
        q4 = (A23 - A32) / (4*q1);
        q2 = (A12 + A21) / (4*q1);
        q3 = (A31 + A13) / (4*q1);
    elseif S3 >= S1 && S3 >= S2 && S3 >= S4
        q2 = 0.5 * sqrt(S3);
        q4 = (A31 - A13) / (4*q2);
        q1 = (A12 + A21) / (4*q2);
        q3 = (A23 + A32) / (4*q2);
    else
        q3 = 0.5 * sqrt(S4);
        q4 = (A12 - A21) / (4*q3);
        q1 = (A31 + A13) / (4*q3);
        q2 = (A23 + A32) / (4*q3);
    end
    
    q_cmd = [q1; q2; q3; q4];
    
    % Normalize and ensure positive scalar
    q_mag = sqrt(q_cmd(1)^2 + q_cmd(2)^2 + q_cmd(3)^2 + q_cmd(4)^2);
    q_cmd = q_cmd / q_mag;
    
    if q_cmd(4) < 0
        q_cmd = -q_cmd;
    end
end


%% ========================================================================
%  8. PD ATTITUDE CONTROLLER
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "PD_Controller"

function tau_cmd = fcn(q, q_cmd, omega, h_wheels, Kp, Kd, tau_max, I, ff_enable)
%#codegen
% PD_CONTROLLER - Quaternion feedback PD controller
%
% Inputs:
%   q        - Current quaternion [4x1]
%   q_cmd    - Desired quaternion [4x1]
%   omega    - Current angular velocity [3x1]
%   h_wheels - Wheel momentum [3x1]
%   Kp       - Proportional gain (scalar)
%   Kd       - Derivative gain (scalar)
%   tau_max  - Max torque (Nm)
%   I        - Inertia matrix [3x3]
%   ff_enable- Feedforward enable (1 or 0)

    % Quaternion error: q_err = q_cmd^(-1) * q
    % q_cmd conjugate
    q_cmd_conj = [-q_cmd(1); -q_cmd(2); -q_cmd(3); q_cmd(4)];
    
    % Quaternion multiplication
    q_err = [q_cmd_conj(4)*q(1) + q_cmd_conj(1)*q(4) + q_cmd_conj(2)*q(3) - q_cmd_conj(3)*q(2);
             q_cmd_conj(4)*q(2) - q_cmd_conj(1)*q(3) + q_cmd_conj(2)*q(4) + q_cmd_conj(3)*q(1);
             q_cmd_conj(4)*q(3) + q_cmd_conj(1)*q(2) - q_cmd_conj(2)*q(1) + q_cmd_conj(3)*q(4);
             q_cmd_conj(4)*q(4) - q_cmd_conj(1)*q(1) - q_cmd_conj(2)*q(2) - q_cmd_conj(3)*q(3)];
    
    % Ensure shortest path
    if q_err(4) < 0
        q_err = -q_err;
    end
    
    % Error quaternion vector part
    q_err_v = [q_err(1); q_err(2); q_err(3)];
    
    % Angular velocity error (assuming omega_cmd = 0)
    omega_err = omega;
    
    % PD control law
    % tau = -Kp * sign(q4) * q_v - Kd * omega_err
    sgn_q4 = 1;
    if q_err(4) < 0
        sgn_q4 = -1;
    end
    
    tau_p = -Kp * sgn_q4 * q_err_v;
    tau_d = -Kd * omega_err;
    
    % Feedforward (gyroscopic compensation)
    if ff_enable > 0.5
        L = I * omega + h_wheels;
        tau_ff = [omega(2)*L(3) - omega(3)*L(2);
                  omega(3)*L(1) - omega(1)*L(3);
                  omega(1)*L(2) - omega(2)*L(1)];
    else
        tau_ff = [0; 0; 0];
    end
    
    tau_cmd = tau_p + tau_d + tau_ff;
    
    % Saturation
    for i = 1:3
        if tau_cmd(i) > tau_max
            tau_cmd(i) = tau_max;
        elseif tau_cmd(i) < -tau_max
            tau_cmd(i) = -tau_max;
        end
    end
end


%% ========================================================================
%  9. ACTUATOR ALLOCATION
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Actuator_Allocation"

function [rw_cmd, mtq_cmd] = fcn(tau_cmd, h_wheels, B_body, tau_max_rw, h_thresh, dump_gain, mtq_max)
%#codegen
% ACTUATOR_ALLOCATION - Distribute torque to reaction wheels and magnetorquers
%
% Inputs:
%   tau_cmd    - Desired control torque [3x1]
%   h_wheels   - Current wheel momentum [3x1]
%   B_body     - Magnetic field in body frame [3x1]
%   tau_max_rw - Max wheel torque (Nm)
%   h_thresh   - Momentum threshold for dumping (Nms)
%   dump_gain  - Momentum dump gain
%   mtq_max    - Max magnetorquer dipole (Am²)

    % Primary: Reaction wheels
    % For orthogonal wheels: rw_cmd_i = -tau_cmd_i
    rw_cmd = -tau_cmd;
    
    % Saturate wheel torque
    for i = 1:3
        if rw_cmd(i) > tau_max_rw
            rw_cmd(i) = tau_max_rw;
        elseif rw_cmd(i) < -tau_max_rw
            rw_cmd(i) = -tau_max_rw;
        end
    end
    
    % Secondary: Magnetorquer for momentum dumping
    mtq_cmd = [0; 0; 0];
    
    B_mag = sqrt(B_body(1)^2 + B_body(2)^2 + B_body(3)^2);
    
    % Check if momentum dumping needed
    h_mag = sqrt(h_wheels(1)^2 + h_wheels(2)^2 + h_wheels(3)^2);
    
    if h_mag > h_thresh && B_mag > 1e-9
        % Momentum dump torque
        tau_dump = -dump_gain * h_wheels;
        
        % Project onto plane perpendicular to B
        dot_tB = tau_dump(1)*B_body(1) + tau_dump(2)*B_body(2) + tau_dump(3)*B_body(3);
        tau_perp = tau_dump - (dot_tB / (B_mag*B_mag)) * B_body;
        
        % Pseudo-inverse: m = (B x tau_perp) / |B|^2
        mtq_cmd = [(B_body(2)*tau_perp(3) - B_body(3)*tau_perp(2)) / (B_mag*B_mag);
                   (B_body(3)*tau_perp(1) - B_body(1)*tau_perp(3)) / (B_mag*B_mag);
                   (B_body(1)*tau_perp(2) - B_body(2)*tau_perp(1)) / (B_mag*B_mag)];
        
        % Saturate
        for i = 1:3
            if mtq_cmd(i) > mtq_max
                mtq_cmd(i) = mtq_max;
            elseif mtq_cmd(i) < -mtq_max
                mtq_cmd(i) = -mtq_max;
            end
        end
    end
end


%% ========================================================================
%  10. WHEEL DYNAMICS
%  ========================================================================
%  Simulink'te: Bu, Integrator bloğu olarak implemente edilir
%  h_wheels_dot = rw_cmd (wheel motor torque = rate of change of momentum)
%  Sadece saturation için MATLAB Function kullanın

function h_sat = fcn(h, h_max)
%#codegen
% WHEEL_SATURATION - Saturate wheel momentum
    
    h_sat = h;
    for i = 1:3
        if h_sat(i) > h_max
            h_sat(i) = h_max;
        elseif h_sat(i) < -h_max
            h_sat(i) = -h_max;
        end
    end
end


%% ========================================================================
%  11. QUATERNION TO EULER ANGLES
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin (visualization için)
%  İsim: "Quat_to_Euler"

function euler_deg = fcn(q)
%#codegen
% QUAT_TO_EULER - Convert quaternion to Euler angles (3-2-1)
%
% Output: euler_deg = [roll; pitch; yaw] in degrees

    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    % Roll (x-axis rotation)
    roll = atan2(2*(q4*q1 + q2*q3), 1 - 2*(q1*q1 + q2*q2));
    
    % Pitch (y-axis rotation)
    sinp = 2*(q4*q2 - q3*q1);
    if sinp > 1
        sinp = 1;
    elseif sinp < -1
        sinp = -1;
    end
    pitch = asin(sinp);
    
    % Yaw (z-axis rotation)
    yaw = atan2(2*(q4*q3 + q1*q2), 1 - 2*(q2*q2 + q3*q3));
    
    % Convert to degrees
    euler_deg = [roll; pitch; yaw] * 57.29577951308232;  % 180/pi
end


%% ========================================================================
%  12. ATTITUDE ERROR COMPUTATION
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin (visualization için)
%  İsim: "Attitude_Error"

function [theta_err_deg, q_err] = fcn(q, q_cmd)
%#codegen
% ATTITUDE_ERROR - Compute attitude error magnitude
%
% Outputs:
%   theta_err_deg - Total attitude error (degrees)
%   q_err         - Error quaternion [4x1]

    % Quaternion error
    q_cmd_conj = [-q_cmd(1); -q_cmd(2); -q_cmd(3); q_cmd(4)];
    
    q_err = [q_cmd_conj(4)*q(1) + q_cmd_conj(1)*q(4) + q_cmd_conj(2)*q(3) - q_cmd_conj(3)*q(2);
             q_cmd_conj(4)*q(2) - q_cmd_conj(1)*q(3) + q_cmd_conj(2)*q(4) + q_cmd_conj(3)*q(1);
             q_cmd_conj(4)*q(3) + q_cmd_conj(1)*q(2) - q_cmd_conj(2)*q(1) + q_cmd_conj(3)*q(4);
             q_cmd_conj(4)*q(4) - q_cmd_conj(1)*q(1) - q_cmd_conj(2)*q(2) - q_cmd_conj(3)*q(3)];
    
    if q_err(4) < 0
        q_err = -q_err;
    end
    
    % Total error angle
    if q_err(4) > 1
        q_err(4) = 1;
    end
    theta_err_deg = 2 * acos(q_err(4)) * 57.29577951308232;
end


%% ========================================================================
%  13. DISTURBANCE TORQUES
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Disturbance_Torques"

function tau_dist = fcn(r_eci, q, I, mu, residual_dipole, B_body)
%#codegen
% DISTURBANCE_TORQUES - Compute environmental disturbance torques
%
% Includes: Gravity gradient, magnetic disturbance

    % DCM
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    A = [q4^2+q1^2-q2^2-q3^2,  2*(q1*q2+q3*q4),      2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4),      q4^2-q1^2+q2^2-q3^2,  2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4),      2*(q2*q3-q1*q4),      q4^2-q1^2-q2^2+q3^2];
    
    % Position magnitude
    r_mag = sqrt(r_eci(1)^2 + r_eci(2)^2 + r_eci(3)^2);
    
    % Nadir direction in body frame
    r_body = A * (r_eci / r_mag);
    
    % Gravity gradient torque (Wertz Eq. 17-31)
    % tau_gg = (3*mu/r^3) * (n x I*n)
    In = I * r_body;
    tau_gg = (3 * mu / r_mag^3) * [r_body(2)*In(3) - r_body(3)*In(2);
                                    r_body(3)*In(1) - r_body(1)*In(3);
                                    r_body(1)*In(2) - r_body(2)*In(1)];
    
    % Magnetic disturbance (residual dipole)
    tau_mag = [residual_dipole(2)*B_body(3) - residual_dipole(3)*B_body(2);
               residual_dipole(3)*B_body(1) - residual_dipole(1)*B_body(3);
               residual_dipole(1)*B_body(2) - residual_dipole(2)*B_body(1)];
    
    tau_dist = tau_gg + tau_mag;
end


%% ========================================================================
%  14. VECTOR ROTATION BY QUATERNION
%  ========================================================================
%  Simulink'te: MATLAB Function bloğu ekleyin
%  İsim: "Rotate_Vector"

function v_body = fcn(v_eci, q)
%#codegen
% ROTATE_VECTOR - Rotate vector from ECI to Body using quaternion
    
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    A = [q4^2+q1^2-q2^2-q3^2,  2*(q1*q2+q3*q4),      2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4),      q4^2-q1^2+q2^2-q3^2,  2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4),      2*(q2*q3-q1*q4),      q4^2-q1^2-q2^2+q3^2];
    
    v_body = A * v_eci;
end
