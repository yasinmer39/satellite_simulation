% Sun Pointing Guidance - Hedef Quaternion Hesabı
% Case Study (g): Aracın yönelim (açı) komutu güneşin konumu üzerinden hesaplanmalıdır
%
% Bu fonksiyon, güneş konumundan hedef quaternion'u hesaplar.
% Amaç: Uzay aracının belirli bir ekseni (örn. +Z) güneşe yönlendirilir.
%
% Referanslar:
%   - Wertz, Chapter 18 & 19 (Attitude Control and Maneuvers)
%   - Markley & Crassidis, "Fundamentals of Spacecraft Attitude Determination"
%   - Wie, "Space Vehicle Dynamics and Control"
%
% Pointing Modları:
%   1. Sun-Pointing: Solar panel normal → güneşe
%   2. Nadir-Pointing + Sun constraint: Z → nadir, X constraint → güneş
%   3. Inertial-Pointing: Sabit inertial yönelim
%
% Hedef Quaternion Hesabı:
%   Body frame'den Reference frame'e dönüşüm quaternion'u
%   q_desired = q_ref_to_body
%
% Author: Mehmet Yasin Meriç
% Date: 2026
% Course: Uzay Sistemleri Vaka Çalışması

function [q_cmd, pointing_info] = sun_pointing_guidance(sun_eci, r_eci, v_eci, mode, params)
% SUN_POINTING_GUIDANCE - Güneş konumundan hedef quaternion hesapla
%
% Girdiler:
%   sun_eci  - Güneş yön vektörü (ECI, unit) [3x1]
%   r_eci    - Uzay aracı pozisyonu (ECI, m) [3x1]
%   v_eci    - Uzay aracı hızı (ECI, m/s) [3x1]
%   mode     - Pointing modu:
%              'sun'    - Pure sun pointing (+Z → sun)
%              'nadir'  - Nadir pointing with sun constraint (+Z → nadir, +X → sun direction)
%              'lvlh'   - LVLH frame alignment
%   params   - Opsiyonel parametreler (struct)
%              .body_axis_to_sun - Güneşe yönlendirilecek body ekseni [3x1]
%              .body_axis_to_nadir - Nadire yönlendirilecek body ekseni [3x1]
%              .offset_angle - Güneş yönelim offset açısı (derece)
%
% Çıktılar:
%   q_cmd        - Hedef quaternion (ECI → Body) [4x1]
%   pointing_info - Bilgi struct'ı (debug için)
%
% Quaternion convention: q = [q1; q2; q3; q4] (scalar-last)

    if nargin < 4 || isempty(mode)
        mode = 'sun';
    end
    
    if nargin < 5 || isempty(params)
        params = struct();
    end
    
    % Default parameters
    if ~isfield(params, 'body_axis_to_sun')
        params.body_axis_to_sun = [0; 0; 1];  % +Z → sun
    end
    if ~isfield(params, 'body_axis_to_nadir')
        params.body_axis_to_nadir = [0; 0; 1];  % +Z → nadir
    end
    if ~isfield(params, 'offset_angle')
        params.offset_angle = 0;  % degrees
    end
    
    % Normalize inputs
    sun_eci = sun_eci / norm(sun_eci);
    r_eci = r_eci(:);
    v_eci = v_eci(:);
    
    % Initialize output
    pointing_info = struct();
    pointing_info.mode = mode;
    pointing_info.sun_eci = sun_eci;
    
    switch lower(mode)
        case 'sun'
            % Pure sun pointing
            % Body axis → Sun direction
            [q_cmd, info] = compute_sun_pointing(sun_eci, params);
            pointing_info.primary_axis = params.body_axis_to_sun;
            
        case 'nadir'
            % Nadir pointing with sun constraint
            % Primary: Body axis → Nadir
            % Secondary: Constrained toward sun
            [q_cmd, info] = compute_nadir_with_sun_constraint(sun_eci, r_eci, v_eci, params);
            pointing_info.primary_axis = params.body_axis_to_nadir;
            pointing_info.nadir_eci = info.nadir_eci;
            
        case 'lvlh'
            % LVLH frame alignment
            [q_cmd, info] = compute_lvlh_frame(r_eci, v_eci);
            pointing_info.lvlh_x = info.x_lvlh;
            pointing_info.lvlh_y = info.y_lvlh;
            pointing_info.lvlh_z = info.z_lvlh;
            
        otherwise
            error('Unknown pointing mode: %s', mode);
    end
    
    % Add common info
    pointing_info.q_cmd = q_cmd;
    
end

%% ==================== PURE SUN POINTING ====================
function [q_cmd, info] = compute_sun_pointing(sun_eci, params)
% Compute quaternion to point body axis toward sun
%
% Wertz Eq. 12-45 benzeri yaklaşım:
% İki vektör çifti kullanarak DCM hesabı

    % Body axis to point toward sun
    body_sun = params.body_axis_to_sun / norm(params.body_axis_to_sun);
    
    % Apply offset if specified
    if params.offset_angle ~= 0
        % Rotate sun vector around an arbitrary perpendicular axis
        perp_axis = cross(sun_eci, [0;0;1]);
        if norm(perp_axis) < 1e-6
            perp_axis = cross(sun_eci, [0;1;0]);
        end
        perp_axis = perp_axis / norm(perp_axis);
        sun_eci = rotate_vector_axis_angle(sun_eci, perp_axis, params.offset_angle * pi/180);
    end
    
    % Reference direction in ECI (target for body axis)
    ref_primary = sun_eci;
    
    % Secondary axis: arbitrary perpendicular (for complete attitude)
    % Use ECI Z or velocity direction as secondary constraint
    ref_secondary = [0; 0; 1];  % ECI Z-axis
    if abs(dot(ref_primary, ref_secondary)) > 0.99
        ref_secondary = [0; 1; 0];  % Use Y if Z is parallel to sun
    end
    
    % Gram-Schmidt to get orthonormal basis in ECI
    z_eci = ref_primary;
    y_eci = cross(z_eci, ref_secondary);
    y_eci = y_eci / norm(y_eci);
    x_eci = cross(y_eci, z_eci);
    
    % Body frame basis (where we want the ECI basis to go)
    z_body = body_sun;
    % Choose arbitrary secondary body axis
    if abs(z_body(3)) < 0.99
        y_body = cross(z_body, [0;0;1]);
    else
        y_body = cross(z_body, [0;1;0]);
    end
    y_body = y_body / norm(y_body);
    x_body = cross(y_body, z_body);
    
    % DCM from ECI to Body
    % A_eci2body * v_eci = v_body
    A_eci = [x_eci, y_eci, z_eci]';  % ECI basis as rows
    A_body = [x_body, y_body, z_body];  % Body basis as columns
    
    A_eci2body = A_body' * A_eci;
    
    % Convert DCM to quaternion
    q_cmd = dcm_to_quat(A_eci2body);
    
    info.A_eci2body = A_eci2body;
    info.sun_angle_offset = params.offset_angle;
end

%% ==================== NADIR POINTING WITH SUN CONSTRAINT ====================
function [q_cmd, info] = compute_nadir_with_sun_constraint(sun_eci, r_eci, v_eci, params)
% Nadir pointing with sun constraint
% Primary: -Z body → nadir (toward Earth center)
% Secondary: +X body constrained toward sun projection

    % Nadir direction (from spacecraft toward Earth center)
    nadir_eci = -r_eci / norm(r_eci);
    
    % Body axis to point toward nadir
    body_nadir = params.body_axis_to_nadir / norm(params.body_axis_to_nadir);
    
    % Sun direction in nadir-perpendicular plane
    sun_proj = sun_eci - dot(sun_eci, nadir_eci) * nadir_eci;
    if norm(sun_proj) < 1e-6
        % Sun is directly overhead/below - use velocity for secondary
        sun_proj = v_eci - dot(v_eci, nadir_eci) * nadir_eci;
    end
    sun_proj = sun_proj / norm(sun_proj);
    
    % Build reference frame in ECI
    % Z_ref → nadir, X_ref → sun projection, Y_ref = Z x X
    z_eci = nadir_eci;
    x_eci = sun_proj;
    y_eci = cross(z_eci, x_eci);
    y_eci = y_eci / norm(y_eci);
    x_eci = cross(y_eci, z_eci);  % Ensure orthogonality
    
    % Build body frame
    % body_nadir axis → nadir
    z_body = body_nadir;
    
    % Choose X body toward sun projection
    % Y body = Z x X
    if abs(z_body(1)) < 0.99
        x_body_init = [1; 0; 0];
    else
        x_body_init = [0; 1; 0];
    end
    y_body = cross(z_body, x_body_init);
    y_body = y_body / norm(y_body);
    x_body = cross(y_body, z_body);
    
    % DCM from ECI to Body
    A_eci = [x_eci, y_eci, z_eci]';
    A_body = [x_body, y_body, z_body];
    
    A_eci2body = A_body' * A_eci;
    
    % Convert DCM to quaternion
    q_cmd = dcm_to_quat(A_eci2body);
    
    info.nadir_eci = nadir_eci;
    info.sun_proj = sun_proj;
    info.A_eci2body = A_eci2body;
end

%% ==================== LVLH FRAME ALIGNMENT ====================
function [q_cmd, info] = compute_lvlh_frame(r_eci, v_eci)
% Compute quaternion for LVLH frame alignment
% LVLH: Local Vertical Local Horizontal
% Z → nadir, X → velocity direction (approximately), Y → orbit normal

    % LVLH basis vectors
    z_lvlh = -r_eci / norm(r_eci);  % Nadir
    h = cross(r_eci, v_eci);        % Angular momentum
    y_lvlh = -h / norm(h);          % Negative orbit normal (for right-hand)
    x_lvlh = cross(y_lvlh, z_lvlh);
    x_lvlh = x_lvlh / norm(x_lvlh);
    
    % DCM from ECI to LVLH
    A_eci2lvlh = [x_lvlh'; y_lvlh'; z_lvlh'];
    
    % Convert DCM to quaternion
    q_cmd = dcm_to_quat(A_eci2lvlh);
    
    info.x_lvlh = x_lvlh;
    info.y_lvlh = y_lvlh;
    info.z_lvlh = z_lvlh;
    info.A_eci2lvlh = A_eci2lvlh;
end

%% ==================== HELPER FUNCTIONS ====================

function q = dcm_to_quat(A)
% Convert DCM to quaternion (scalar-last convention)
% Shepperd's method for numerical stability
%
% Reference: Wertz Eq. 12-13

    tr = trace(A);
    
    % Shepperd's method
    S = [1 + tr;
         1 + 2*A(1,1) - tr;
         1 + 2*A(2,2) - tr;
         1 + 2*A(3,3) - tr];
    
    [~, idx] = max(S);
    
    switch idx
        case 1
            q4 = 0.5 * sqrt(S(1));
            q1 = (A(2,3) - A(3,2)) / (4*q4);
            q2 = (A(3,1) - A(1,3)) / (4*q4);
            q3 = (A(1,2) - A(2,1)) / (4*q4);
        case 2
            q1 = 0.5 * sqrt(S(2));
            q4 = (A(2,3) - A(3,2)) / (4*q1);
            q2 = (A(1,2) + A(2,1)) / (4*q1);
            q3 = (A(3,1) + A(1,3)) / (4*q1);
        case 3
            q2 = 0.5 * sqrt(S(3));
            q4 = (A(3,1) - A(1,3)) / (4*q2);
            q1 = (A(1,2) + A(2,1)) / (4*q2);
            q3 = (A(2,3) + A(3,2)) / (4*q2);
        case 4
            q3 = 0.5 * sqrt(S(4));
            q4 = (A(1,2) - A(2,1)) / (4*q3);
            q1 = (A(3,1) + A(1,3)) / (4*q3);
            q2 = (A(2,3) + A(3,2)) / (4*q3);
    end
    
    q = [q1; q2; q3; q4];
    
    % Ensure scalar component is positive (avoid unwinding)
    if q(4) < 0
        q = -q;
    end
    
    % Normalize
    q = q / norm(q);
end

function v_rot = rotate_vector_axis_angle(v, axis, angle)
% Rotate vector v around axis by angle (Rodrigues' formula)
    axis = axis / norm(axis);
    v_rot = v * cos(angle) + cross(axis, v) * sin(angle) + ...
            axis * dot(axis, v) * (1 - cos(angle));
end
