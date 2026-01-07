% Quaternion Utility Functions
% Attitude determination için gerekli quaternion işlemleri.
%
% Quaternion Convention: q = [q1, q2, q3, q4] where q4 is scalar part
%   q = [e*sin(θ/2); cos(θ/2)]
%   e: rotation axis (unit vector)
%   θ: rotation angle
%
% Referanslar:
%   - Wertz Chapter 12.1 "Parameterization of the Attitude"
%   - Shuster (1993), "A Survey of Attitude Representations"
%   - Markley & Crassidis (2014), "Fundamentals of Spacecraft Attitude Determination"

%% ==================== QUATERNION OPERATIONS ====================

function r = quat_multiply(p, q)
% Quaternion multiplication: r = p ⊗ q
% Represents rotation p followed by rotation q
%
% p, q: quaternions [4x1], scalar-last convention
% r: result quaternion [4x1]

    p1 = p(1); p2 = p(2); p3 = p(3); p4 = p(4);
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    r = [p4*q1 + p1*q4 + p2*q3 - p3*q2;
         p4*q2 - p1*q3 + p2*q4 + p3*q1;
         p4*q3 + p1*q2 - p2*q1 + p3*q4;
         p4*q4 - p1*q1 - p2*q2 - p3*q3];
end


function q_conj = quat_conjugate(q)
% Quaternion conjugate (inverse for unit quaternions)
% q_conj = [-q_vec; q_scalar]
%
% q: quaternion [4x1]
% q_conj: conjugate quaternion [4x1]

    q_conj = [-q(1); -q(2); -q(3); q(4)];
end


function q_inv = quat_inverse(q)
% Quaternion inverse
% For unit quaternions, inverse = conjugate
%
% q: quaternion [4x1]
% q_inv: inverse quaternion [4x1]

    q_norm_sq = sum(q.^2);
    q_inv = quat_conjugate(q) / q_norm_sq;
end


function q_norm = quat_normalize(q)
% Quaternion normalization
%
% q: quaternion [4x1]
% q_norm: unit quaternion [4x1]

    q_norm = q / norm(q);
end


function v_rot = quat_rotate(q, v)
% Rotate vector v by quaternion q
% v_rot = q ⊗ v ⊗ q*
%
% q: quaternion [4x1]
% v: vector [3x1]
% v_rot: rotated vector [3x1]

    v = v(:);
    
    % Expand v to quaternion form [v; 0]
    v_quat = [v; 0];
    
    % v_rot = q * v * q_conj
    q_conj = quat_conjugate(q);
    temp = quat_multiply(q, v_quat);
    v_rot_quat = quat_multiply(temp, q_conj);
    
    v_rot = v_rot_quat(1:3);
end


%% ==================== CONVERSIONS ====================

function A = quat_to_dcm(q)
% Quaternion to Direction Cosine Matrix (DCM)
% Wertz Eq. 12-13a
%
% q: quaternion [4x1], scalar-last
% A: DCM [3x3], body <- reference transformation

    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    A = [q1^2-q2^2-q3^2+q4^2,   2*(q1*q2+q3*q4),       2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4),       -q1^2+q2^2-q3^2+q4^2,  2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4),       2*(q2*q3-q1*q4),       -q1^2-q2^2+q3^2+q4^2];
end


function q = dcm_to_quat(A)
% DCM to Quaternion conversion
% Shuster method for numerical stability
%
% A: DCM [3x3]
% q: quaternion [4x1], scalar-last

    tr = trace(A);
    
    % K matrix (Shuster method)
    K = [A(1,1)-A(2,2)-A(3,3),  A(2,1)+A(1,2),          A(3,1)+A(1,3),          A(2,3)-A(3,2);
         A(2,1)+A(1,2),          A(2,2)-A(1,1)-A(3,3),  A(3,2)+A(2,3),          A(3,1)-A(1,3);
         A(3,1)+A(1,3),          A(3,2)+A(2,3),          A(3,3)-A(1,1)-A(2,2),  A(1,2)-A(2,1);
         A(2,3)-A(3,2),          A(3,1)-A(1,3),          A(1,2)-A(2,1),          tr] / 3;
    
    [V, D] = eig(K);
    [~, idx] = max(diag(D));
    q = V(:, idx);
    q = q / norm(q);
    
    if q(4) < 0
        q = -q;
    end
end


function euler = quat_to_euler(q)
% Quaternion to Euler angles (3-2-1 sequence: yaw-pitch-roll)
%
% q: quaternion [4x1], scalar-last
% euler: [roll, pitch, yaw] in radians [3x1]

    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    
    % Roll (x-axis rotation)
    roll = atan2(2*(q4*q1 + q2*q3), 1 - 2*(q1^2 + q2^2));
    
    % Pitch (y-axis rotation)
    sinp = 2*(q4*q2 - q3*q1);
    if abs(sinp) >= 1
        pitch = sign(sinp) * pi/2;  % Gimbal lock
    else
        pitch = asin(sinp);
    end
    
    % Yaw (z-axis rotation)
    yaw = atan2(2*(q4*q3 + q1*q2), 1 - 2*(q2^2 + q3^2));
    
    euler = [roll; pitch; yaw];
end


function q = euler_to_quat(euler)
% Euler angles to Quaternion (3-2-1 sequence)
%
% euler: [roll, pitch, yaw] in radians [3x1]
% q: quaternion [4x1], scalar-last

    roll = euler(1);
    pitch = euler(2);
    yaw = euler(3);
    
    cr = cos(roll/2);  sr = sin(roll/2);
    cp = cos(pitch/2); sp = sin(pitch/2);
    cy = cos(yaw/2);   sy = sin(yaw/2);
    
    q = [sr*cp*cy - cr*sp*sy;
         cr*sp*cy + sr*cp*sy;
         cr*cp*sy - sr*sp*cy;
         cr*cp*cy + sr*sp*sy];
    
    q = q / norm(q);
end


function [axis, angle] = quat_to_axis_angle(q)
% Quaternion to Axis-Angle representation
% Wertz Eq. 12-7
%
% q: quaternion [4x1], scalar-last
% axis: rotation axis (unit vector) [3x1]
% angle: rotation angle (rad)

    q = q / norm(q);
    
    angle = 2 * acos(q(4));
    
    if abs(angle) < 1e-10
        axis = [0; 0; 1];  % Arbitrary axis for zero rotation
    else
        axis = q(1:3) / sin(angle/2);
        axis = axis / norm(axis);
    end
end


function q = axis_angle_to_quat(axis, angle)
% Axis-Angle to Quaternion conversion
% Wertz Eq. 12-11
%
% axis: rotation axis (unit vector) [3x1]
% angle: rotation angle (rad)
% q: quaternion [4x1], scalar-last

    axis = axis(:) / norm(axis);
    
    q = [axis * sin(angle/2);
         cos(angle/2)];
end


%% ==================== ATTITUDE ERROR ====================

function delta_q = quat_error(q_est, q_true)
% Quaternion error (multiplicative)
% delta_q = q_true ⊗ q_est^(-1)
%
% q_est: estimated quaternion [4x1]
% q_true: true quaternion [4x1]
% delta_q: error quaternion [4x1]

    delta_q = quat_multiply(q_true, quat_inverse(q_est));
    
    % Ensure scalar part is positive
    if delta_q(4) < 0
        delta_q = -delta_q;
    end
end


function angle = quat_angle_error(q_est, q_true)
% Angular error between two quaternions (deg)
%
% q_est: estimated quaternion [4x1]
% q_true: true quaternion [4x1]
% angle: angular error (deg)

    delta_q = quat_error(q_est, q_true);
    
    % Rotation angle
    angle = 2 * acos(abs(delta_q(4)));
    angle = rad2deg(angle);
end


%% ==================== INTERPOLATION ====================

function q_interp = quat_slerp(q1, q2, t)
% Spherical Linear Interpolation (SLERP)
%
% q1, q2: quaternions [4x1]
% t: interpolation parameter [0, 1]
% q_interp: interpolated quaternion [4x1]

    q1 = q1 / norm(q1);
    q2 = q2 / norm(q2);
    
    dot_prod = sum(q1 .* q2);
    
    % Ensure shortest path
    if dot_prod < 0
        q2 = -q2;
        dot_prod = -dot_prod;
    end
    
    % Threshold for linear interpolation
    if dot_prod > 0.9995
        q_interp = q1 + t * (q2 - q1);
        q_interp = q_interp / norm(q_interp);
        return;
    end
    
    theta_0 = acos(dot_prod);
    theta = theta_0 * t;
    
    q_interp = (sin(theta_0 - theta) / sin(theta_0)) * q1 + ...
               (sin(theta) / sin(theta_0)) * q2;
    q_interp = q_interp / norm(q_interp);
end
