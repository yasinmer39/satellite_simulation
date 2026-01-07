% Actuator Allocation - Kontrol Torku Dağıtımı
% Case Study (g): Tasarladığımız eyleyicilere çıktı göndermeli
%
% Bu fonksiyon, kontrolcüden gelen toplam tork komutunu
% reaction wheel ve magnetorquer eyleyicilerine dağıtır.
%
% Eyleyiciler:
%   1. Reaction Wheels (3 adet, orthogonal) - Hızlı, hassas kontrol
%   2. Magnetorquers (3 adet, orthogonal) - Momentum dumping, düşük tork
%
% Allocation Stratejisi:
%   - High frequency torque → Reaction wheels
%   - DC/secular torque → Magnetorquers (momentum dumping)
%   - Saturation handling → Redistribution
%
% Referanslar:
%   - Wertz Chapter 18.2 (Reaction Wheel Systems)
%   - Sidi, "Spacecraft Dynamics and Control", Chapter 7
%
% Author: Mehmet Yasin Meriç
% Date: 2026
% Course: Uzay Sistemleri Vaka Çalışması

function [wheel_cmd, mtq_cmd, alloc_info] = actuator_allocation(tau_cmd, h_wheels, B_body, params)
% ACTUATOR_ALLOCATION - Kontrol torkunu eyleyicilere dağıt
%
% Girdiler:
%   tau_cmd  - İstenen kontrol torku (body frame, Nm) [3x1]
%   h_wheels - Mevcut wheel angular momentum (body frame, Nms) [3x1]
%   B_body   - Manyetik alan vektörü (body frame, T) [3x1]
%   params   - Eyleyici parametreleri (struct):
%              .wheel_max_torque - Max wheel torque (Nm) [3x1] veya scalar
%              .wheel_max_momentum - Max wheel momentum (Nms) [3x1] veya scalar
%              .mtq_max_dipole - Max magnetorquer dipole (Am²) [3x1] veya scalar
%              .momentum_dump_gain - Momentum dumping gain
%              .momentum_threshold - Threshold for momentum dumping
%              .allocation_mode - 'wheel_only', 'mtq_only', 'hybrid'
%
% Çıktılar:
%   wheel_cmd  - Wheel torque commands (Nm) [3x1]
%   mtq_cmd    - Magnetorquer dipole commands (Am²) [3x1]
%   alloc_info - Allocation bilgileri (struct)
%
% Not: wheel_cmd, wheel motor'a gönderilecek torque command'dır.
%      Pozitif değer wheel'ı hızlandırır, spacecraft'a negatif torque uygular.
%      tau_applied = -wheel_cmd (reaction)

    %% Input validation and defaults
    tau_cmd = tau_cmd(:);
    h_wheels = h_wheels(:);
    B_body = B_body(:);
    
    % Default parameters
    if ~isfield(params, 'wheel_max_torque')
        params.wheel_max_torque = 0.01;  % 10 mNm (small satellite RW)
    end
    if ~isfield(params, 'wheel_max_momentum')
        params.wheel_max_momentum = 0.4;  % 0.4 Nms
    end
    if ~isfield(params, 'mtq_max_dipole')
        params.mtq_max_dipole = 5;  % 5 Am²
    end
    if ~isfield(params, 'momentum_dump_gain')
        params.momentum_dump_gain = 0.1;
    end
    if ~isfield(params, 'momentum_threshold')
        params.momentum_threshold = 0.3;  % 75% of max
    end
    if ~isfield(params, 'allocation_mode')
        params.allocation_mode = 'hybrid';
    end
    
    % Convert scalar limits to vectors
    if isscalar(params.wheel_max_torque)
        wheel_max_tau = params.wheel_max_torque * ones(3,1);
    else
        wheel_max_tau = params.wheel_max_torque(:);
    end
    
    if isscalar(params.wheel_max_momentum)
        wheel_max_h = params.wheel_max_momentum * ones(3,1);
    else
        wheel_max_h = params.wheel_max_momentum(:);
    end
    
    if isscalar(params.mtq_max_dipole)
        mtq_max_m = params.mtq_max_dipole * ones(3,1);
    else
        mtq_max_m = params.mtq_max_dipole(:);
    end
    
    %% Initialize outputs
    wheel_cmd = zeros(3, 1);
    mtq_cmd = zeros(3, 1);
    
    alloc_info.mode = params.allocation_mode;
    alloc_info.wheel_saturation = false;
    alloc_info.mtq_active = false;
    alloc_info.momentum_dump = false;
    
    %% Allocation based on mode
    switch lower(params.allocation_mode)
        case 'wheel_only'
            % All torque to reaction wheels
            wheel_cmd = allocate_to_wheels(tau_cmd, wheel_max_tau, h_wheels, wheel_max_h);
            
        case 'mtq_only'
            % All torque to magnetorquers
            mtq_cmd = allocate_to_magnetorquers(tau_cmd, B_body, mtq_max_m);
            alloc_info.mtq_active = true;
            
        case 'hybrid'
            % Primary: Reaction wheels
            % Secondary: Magnetorquers for momentum dumping
            
            % Step 1: Allocate primary torque to wheels
            wheel_cmd = allocate_to_wheels(tau_cmd, wheel_max_tau, h_wheels, wheel_max_h);
            
            % Step 2: Check if momentum dumping is needed
            h_norm = abs(h_wheels);
            need_dump = any(h_norm > params.momentum_threshold * wheel_max_h);
            
            if need_dump && norm(B_body) > 1e-9
                % Momentum dumping torque
                % τ_dump = -K * h_wheels (to reduce stored momentum)
                tau_dump = -params.momentum_dump_gain * h_wheels;
                
                % Allocate dump torque to magnetorquers
                mtq_cmd = allocate_to_magnetorquers(tau_dump, B_body, mtq_max_m);
                
                alloc_info.momentum_dump = true;
                alloc_info.mtq_active = true;
                alloc_info.tau_dump = tau_dump;
            end
            
            % Step 3: Handle wheel torque saturation
            tau_residual = tau_cmd - (-wheel_cmd);  % Unmet torque request
            if norm(tau_residual) > 1e-9 && norm(B_body) > 1e-9
                % Use magnetorquers to supplement
                mtq_cmd_supp = allocate_to_magnetorquers(tau_residual, B_body, mtq_max_m);
                mtq_cmd = mtq_cmd + mtq_cmd_supp;
                alloc_info.mtq_active = true;
            end
            
        otherwise
            error('Unknown allocation mode: %s', params.allocation_mode);
    end
    
    %% Saturate magnetorquer commands
    for i = 1:3
        if abs(mtq_cmd(i)) > mtq_max_m(i)
            mtq_cmd(i) = sign(mtq_cmd(i)) * mtq_max_m(i);
        end
    end
    
    %% Compute actual applied torques
    tau_wheel = -wheel_cmd;  % Reaction torque
    tau_mtq = cross(mtq_cmd, B_body);  % τ = m × B
    
    alloc_info.tau_wheel = tau_wheel;
    alloc_info.tau_mtq = tau_mtq;
    alloc_info.tau_total = tau_wheel + tau_mtq;
    alloc_info.tau_requested = tau_cmd;
    alloc_info.tau_error = tau_cmd - alloc_info.tau_total;
    
    % Check saturation
    alloc_info.wheel_saturation = any(abs(wheel_cmd) >= wheel_max_tau - 1e-10);
    alloc_info.h_wheels_pct = abs(h_wheels) ./ wheel_max_h * 100;
    
end

%% ==================== WHEEL ALLOCATION ====================
function wheel_cmd = allocate_to_wheels(tau_cmd, tau_max, h_current, h_max)
% Allocate torque to orthogonal reaction wheels
%
% For orthogonal 3-wheel configuration:
% wheel_cmd_i = -tau_cmd_i (direct mapping)
%
% Note: wheel_cmd is the torque applied TO the wheel
%       This creates reaction torque ON the spacecraft = -wheel_cmd

    wheel_cmd = -tau_cmd;  % Direct allocation for orthogonal wheels
    
    % Check for momentum saturation risk
    % If wheel is near saturation and command would push it further, limit
    for i = 1:3
        % Predicted momentum after applying torque for some time
        % This is a simple check - more sophisticated would use time horizon
        
        % Saturate torque command
        if abs(wheel_cmd(i)) > tau_max(i)
            wheel_cmd(i) = sign(wheel_cmd(i)) * tau_max(i);
        end
        
        % Additional check: don't command torque that would exceed momentum limit
        % If momentum is near max and command would increase it, limit
        if abs(h_current(i)) > 0.9 * h_max(i)
            if sign(wheel_cmd(i)) == sign(h_current(i))
                % Would increase momentum magnitude - limit
                wheel_cmd(i) = wheel_cmd(i) * 0.5;
            end
        end
    end
end

%% ==================== MAGNETORQUER ALLOCATION ====================
function mtq_cmd = allocate_to_magnetorquers(tau_cmd, B_body, m_max)
% Allocate torque to magnetorquers
%
% Magnetorquer physics:
%   τ = m × B
%
% Given desired τ and measured B, find m:
%   This is an underdetermined problem (2 DOF from 3 unknowns)
%   We use pseudo-inverse: m = B× τ / |B|²
%   But τ component parallel to B cannot be generated!
%
% Reference: Wertz Section 19.3 (Magnetic Control)

    B_norm = norm(B_body);
    
    if B_norm < 1e-9
        mtq_cmd = zeros(3, 1);
        return;
    end
    
    % Only the component of τ perpendicular to B can be generated
    tau_perp = tau_cmd - dot(tau_cmd, B_body) / B_norm^2 * B_body;
    
    % Pseudo-inverse solution
    % m = (B × τ) / |B|²
    mtq_cmd = cross(B_body, tau_perp) / B_norm^2;
    
    % Saturate
    for i = 1:3
        if abs(mtq_cmd(i)) > m_max(i)
            mtq_cmd(i) = sign(mtq_cmd(i)) * m_max(i);
        end
    end
end

%% ==================== WHEEL DYNAMICS ====================
function [h_new, omega_wheel_new] = wheel_dynamics(h_current, wheel_cmd, dt, params)
% WHEEL_DYNAMICS - Reaction wheel momentum update
%
% ḣ_wheel = τ_motor (wheel_cmd)
% τ_applied_to_sc = -τ_motor = -wheel_cmd
%
% Girdiler:
%   h_current   - Current wheel momentum [3x1]
%   wheel_cmd   - Motor torque command [3x1]
%   dt          - Time step
%   params      - Wheel parameters
%
% Çıktılar:
%   h_new       - Updated wheel momentum [3x1]
%   omega_wheel_new - Wheel speeds (rad/s) [3x1]

    if ~isfield(params, 'I_wheel')
        params.I_wheel = 1e-4;  % Wheel inertia (kg·m²)
    end
    if ~isfield(params, 'wheel_max_momentum')
        params.wheel_max_momentum = 0.4;
    end
    
    if isscalar(params.I_wheel)
        I_wheel = params.I_wheel * ones(3,1);
    else
        I_wheel = params.I_wheel(:);
    end
    
    if isscalar(params.wheel_max_momentum)
        h_max = params.wheel_max_momentum * ones(3,1);
    else
        h_max = params.wheel_max_momentum(:);
    end
    
    % Integrate wheel momentum
    h_new = h_current + wheel_cmd * dt;
    
    % Saturate momentum
    for i = 1:3
        if abs(h_new(i)) > h_max(i)
            h_new(i) = sign(h_new(i)) * h_max(i);
        end
    end
    
    % Compute wheel speeds
    omega_wheel_new = h_new ./ I_wheel;
end
