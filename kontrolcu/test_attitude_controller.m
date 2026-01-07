% Attitude Controller Test Script
% Case Study (g): PID kontrolcü tasarlayınız
%
% Bu script, sun pointing kontrolcüsünü test eder:
%   1. Sun pointing guidance
%   2. PD attitude controller
%   3. Actuator allocation
%   4. Closed-loop simulation
%   5. Performance analysis
%
% Author: Mehmet Yasin Meriç
% Date: 2026
% Course: Uzay Sistemleri Vaka Çalışması

clear; clc; close all;

fprintf('=================================================================\n');
fprintf('   CASE STUDY (g): ATTITUDE CONTROLLER TEST\n');
fprintf('   60 kg LEO Satellite - Sun Pointing Control\n');
fprintf('=================================================================\n\n');

%% ==================== SPACECRAFT PARAMETERS ====================
fprintf('1. Spacecraft Parametreleri\n');
fprintf('-----------------------------------------------------------------\n');

% Mass properties
spacecraft.mass = 60;  % kg
spacecraft.I = diag([5.4, 5.4, 4.32]);  % kg·m² (from earlier design)

% Orbit
orbit.altitude = 1000e3;  % m
orbit.R_earth = 6378137;  % m
orbit.mu = 3.986004418e14;  % m³/s²
orbit.a = orbit.R_earth + orbit.altitude;
orbit.n = sqrt(orbit.mu / orbit.a^3);  % rad/s
orbit.period = 2*pi / orbit.n;
orbit.inclination = 60 * pi/180;  % rad

fprintf('   - Mass: %.1f kg\n', spacecraft.mass);
fprintf('   - Inertia: [%.2f, %.2f, %.2f] kg·m²\n', spacecraft.I(1,1), spacecraft.I(2,2), spacecraft.I(3,3));
fprintf('   - Orbit altitude: %.0f km\n', orbit.altitude/1000);
fprintf('   - Orbit period: %.1f min\n', orbit.period/60);

%% ==================== ACTUATOR PARAMETERS ====================
fprintf('\n2. Eyleyici Parametreleri\n');
fprintf('-----------------------------------------------------------------\n');

% Reaction Wheels (based on earlier design - RW400 style)
actuators.wheel_max_torque = 0.01;      % 10 mNm
actuators.wheel_max_momentum = 0.4;      % 0.4 Nms
actuators.wheel_inertia = 1e-4;         % kg·m²
actuators.wheel_max_speed = 6000 * 2*pi/60;  % rad/s

% Magnetorquers (based on earlier design)
actuators.mtq_max_dipole = 5;            % 5 Am²
actuators.mtq_residual = 0.5;            % 0.5 Am² residual

fprintf('   Reaction Wheels (3x orthogonal):\n');
fprintf('      - Max torque: %.1f mNm\n', actuators.wheel_max_torque*1000);
fprintf('      - Max momentum: %.2f Nms\n', actuators.wheel_max_momentum);
fprintf('   Magnetorquers (3x orthogonal):\n');
fprintf('      - Max dipole: %.1f Am²\n', actuators.mtq_max_dipole);

%% ==================== CONTROLLER DESIGN ====================
fprintf('\n3. Kontrolcü Tasarımı (PD)\n');
fprintf('-----------------------------------------------------------------\n');

% Design specifications
% NOT: Büyük açı manevraları için daha agresif tasarım
ctrl_spec.settling_time = 60;      % seconds (hedef - daha hızlı)
ctrl_spec.overshoot = 10;          % percent
ctrl_spec.damping_ratio = 0.707;   % ζ = 1/√2 - optimal damping

% Natural frequency from settling time
% For 2% settling: ts ≈ 4/(ζ*ωn)
ctrl_spec.omega_n = 4 / (ctrl_spec.damping_ratio * ctrl_spec.settling_time);

% Compute PD gains
% Wertz Eq. 18-12:
%   Kp = I * ωn²
%   Kd = 2 * ζ * ωn * I
I_avg = mean(diag(spacecraft.I));

% Torque-limited gain design
% Max quaternion error component ≈ 0.5 for large angles
% tau_max = Kp * q_err_max → Kp = tau_max / q_err_max
q_err_max = 0.4;  % Approximate max quaternion error component
ctrl_params.Kp = actuators.wheel_max_torque / q_err_max * 0.9;  % 90% of limit

% Kd from damping ratio: Kd = 2*ζ*sqrt(Kp*I)
ctrl_params.Kd = 2 * ctrl_spec.damping_ratio * sqrt(ctrl_params.Kp * I_avg);

ctrl_params.I = spacecraft.I;
ctrl_params.tau_max = actuators.wheel_max_torque;
ctrl_params.feedforward = true;

% Recalculate expected performance
omega_n_actual = sqrt(ctrl_params.Kp / I_avg);
zeta_actual = ctrl_params.Kd / (2 * sqrt(ctrl_params.Kp * I_avg));
ts_expected = 4 / (zeta_actual * omega_n_actual);

fprintf('   Design Specifications:\n');
fprintf('      - Target settling time: %.0f s\n', ctrl_spec.settling_time);
fprintf('      - Target damping ratio: %.3f\n', ctrl_spec.damping_ratio);
fprintf('   Computed Gains (torque-limited design):\n');
fprintf('      - Kp = %.4f Nm/rad\n', ctrl_params.Kp);
fprintf('      - Kd = %.4f Nm/(rad/s)\n', ctrl_params.Kd);
fprintf('   Actual closed-loop parameters:\n');
fprintf('      - ωn = %.4f rad/s\n', omega_n_actual);
fprintf('      - ζ  = %.3f\n', zeta_actual);
fprintf('      - Expected ts (linear) = %.1f s\n', ts_expected);
fprintf('   NOTE: Large angle maneuvers will be slower due to saturation\n');

%% ==================== ALLOCATION PARAMETERS ====================
fprintf('\n4. Actuator Allocation Parametreleri\n');
fprintf('-----------------------------------------------------------------\n');

alloc_params.wheel_max_torque = actuators.wheel_max_torque;
alloc_params.wheel_max_momentum = actuators.wheel_max_momentum;
alloc_params.mtq_max_dipole = actuators.mtq_max_dipole;
alloc_params.momentum_dump_gain = 0.05;
alloc_params.momentum_threshold = 0.7 * actuators.wheel_max_momentum;
alloc_params.allocation_mode = 'hybrid';

fprintf('   Mode: %s\n', alloc_params.allocation_mode);
fprintf('   Momentum dump threshold: %.0f%%\n', alloc_params.momentum_threshold/actuators.wheel_max_momentum*100);

%% ==================== SIMULATION SETUP ====================
fprintf('\n5. Simülasyon\n');
fprintf('-----------------------------------------------------------------\n');

% Simulation parameters
t_sim = 600;        % 10 minutes
dt = 0.1;           % 10 Hz
t = 0:dt:t_sim;
N = length(t);

% Initial conditions
% Test Case 1: Moderate initial error (more realistic)
q_init = euler_to_quat([15; -10; 20] * pi/180);  % ~30° total initial error
omega_init = [0.2; -0.1; 0.15] * pi/180;         % Small initial rate (deg/s)
h_wheels_init = [0; 0; 0];                        % Zero wheel momentum

% % Test Case 2: Large initial error (stress test)
% q_init = euler_to_quat([25; -15; 40] * pi/180);  % ~50° total initial error
% omega_init = [0.5; -0.3; 0.2] * pi/180;

% Storage
q_history = zeros(4, N);
omega_history = zeros(3, N);
q_cmd_history = zeros(4, N);
tau_cmd_history = zeros(3, N);
wheel_cmd_history = zeros(3, N);
mtq_cmd_history = zeros(3, N);
h_wheels_history = zeros(3, N);
sun_eci_history = zeros(3, N);
theta_err_history = zeros(1, N);

% Initial state
q_history(:,1) = q_init;
omega_history(:,1) = omega_init;
h_wheels_history(:,1) = h_wheels_init;

% UTC start time (summer solstice for good sun)
utc_start = datetime(2026, 6, 21, 12, 0, 0);

fprintf('   Duration: %.0f s (%.1f min)\n', t_sim, t_sim/60);
fprintf('   Initial attitude error: [25, -15, 40] deg\n');
fprintf('   Starting simulation...\n');

%% ==================== MAIN SIMULATION LOOP ====================
tic;

for k = 1:N
    % Current time
    t_now = t(k);
    utc_now = utc_start + seconds(t_now);
    
    % ===== ENVIRONMENT =====
    % Orbital position (circular orbit)
    M = orbit.n * t_now;  % Mean anomaly
    r_orbital = orbit.a * [cos(M); sin(M); 0];
    R_inc = [1, 0, 0; 0, cos(orbit.inclination), -sin(orbit.inclination); ...
             0, sin(orbit.inclination), cos(orbit.inclination)];
    r_eci = R_inc * r_orbital;
    v_eci = R_inc * orbit.a * orbit.n * [-sin(M); cos(M); 0];
    
    % Sun position
    sun_eci = compute_sun_position(utc_now);
    sun_eci_history(:,k) = sun_eci;
    
    % Magnetic field (dipole model)
    B_eci = compute_mag_field(r_eci);
    A_current = quat_to_dcm(q_history(:,k));
    B_body = A_current * B_eci;
    
    % ===== GUIDANCE =====
    % Compute desired quaternion (sun pointing)
    guidance_params.body_axis_to_sun = [0; 0; 1];  % +Z to sun
    [q_cmd, ~] = sun_pointing_guidance(sun_eci, r_eci, v_eci, 'sun', guidance_params);
    q_cmd_history(:,k) = q_cmd;
    
    % Desired angular velocity (nadir tracking would have non-zero)
    omega_cmd = [0; 0; 0];  % Inertial pointing
    
    % ===== CONTROLLER =====
    ctrl_params.h_wheels = h_wheels_history(:,max(1,k-1));
    [tau_cmd, ctrl_info] = pd_attitude_controller(...
        q_history(:,k), q_cmd, omega_history(:,k), omega_cmd, ctrl_params);
    tau_cmd_history(:,k) = tau_cmd;
    theta_err_history(k) = ctrl_info.theta_error_deg;
    
    % ===== ALLOCATION =====
    [wheel_cmd, mtq_cmd, ~] = actuator_allocation(...
        tau_cmd, h_wheels_history(:,max(1,k-1)), B_body, alloc_params);
    wheel_cmd_history(:,k) = wheel_cmd;
    mtq_cmd_history(:,k) = mtq_cmd;
    
    % ===== DYNAMICS =====
    if k < N
        % Actual applied torques
        tau_wheel = -wheel_cmd;  % Reaction
        tau_mtq = cross(mtq_cmd, B_body);
        tau_total = tau_wheel + tau_mtq;
        
        % Add disturbance torques (simplified)
        tau_dist = [1e-7; -0.5e-7; 0.8e-7];  % Small constant disturbance
        
        % Euler's equation: I*ω̇ + ω × (I*ω + h_w) = τ
        h_wheels = h_wheels_history(:,k);
        omega = omega_history(:,k);
        L = spacecraft.I * omega + h_wheels;
        omega_dot = spacecraft.I \ (tau_total + tau_dist - cross(omega, L));
        
        % Integrate angular velocity
        omega_history(:,k+1) = omega + omega_dot * dt;
        
        % Quaternion kinematics: q̇ = 0.5 * Ω(ω) * q
        q_dot = 0.5 * quat_omega_matrix(omega) * q_history(:,k);
        q_new = q_history(:,k) + q_dot * dt;
        q_history(:,k+1) = q_new / norm(q_new);
        
        % Wheel momentum update
        h_wheels_history(:,k+1) = h_wheels + wheel_cmd * dt;
        
        % Saturate wheel momentum
        for i = 1:3
            if abs(h_wheels_history(i,k+1)) > actuators.wheel_max_momentum
                h_wheels_history(i,k+1) = sign(h_wheels_history(i,k+1)) * actuators.wheel_max_momentum;
            end
        end
    end
end

sim_time = toc;
fprintf('   Simulation completed in %.2f s\n', sim_time);

%% ==================== RESULTS ANALYSIS ====================
fprintf('\n6. Sonuçlar\n');
fprintf('-----------------------------------------------------------------\n');

% Convert to Euler angles for plotting
euler_history = zeros(3, N);
euler_cmd_history = zeros(3, N);
for k = 1:N
    euler_history(:,k) = quat_to_euler(q_history(:,k)) * 180/pi;
    euler_cmd_history(:,k) = quat_to_euler(q_cmd_history(:,k)) * 180/pi;
end

% Performance metrics
settling_idx = find(theta_err_history < 1, 1, 'first');  % 1 deg threshold
if ~isempty(settling_idx)
    settling_time_actual = t(settling_idx);
else
    settling_time_actual = t_sim;
end

ss_error_mean = mean(theta_err_history(end-100:end));
ss_error_std = std(theta_err_history(end-100:end));
max_overshoot = max(theta_err_history);
max_torque_used = max(abs(tau_cmd_history), [], 2);
max_momentum = max(abs(h_wheels_history), [], 2);

fprintf('   Performance Metrics:\n');
fprintf('      - Settling time (1°): %.1f s\n', settling_time_actual);
fprintf('      - Steady-state error: %.3f ± %.3f deg\n', ss_error_mean, ss_error_std);
fprintf('      - Peak attitude error: %.2f deg\n', max_overshoot);
fprintf('      - Max torque used: [%.4f, %.4f, %.4f] Nm\n', max_torque_used);
fprintf('      - Max wheel momentum: [%.3f, %.3f, %.3f] Nms\n', max_momentum);
fprintf('      - Momentum utilization: %.0f%%\n', max(max_momentum)/actuators.wheel_max_momentum*100);

%% ==================== PLOTS ====================
fprintf('\n7. Grafikler\n');
fprintf('-----------------------------------------------------------------\n');

% Figure 1: Attitude Response
figure('Name', 'Attitude Control Response', 'Position', [50, 50, 1400, 900]);

subplot(3,2,1);
plot(t, euler_history(1,:), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, euler_cmd_history(1,:), 'r--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Roll (deg)');
title('Roll Angle');
legend('Actual', 'Command');
grid on;

subplot(3,2,3);
plot(t, euler_history(2,:), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, euler_cmd_history(2,:), 'r--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Pitch (deg)');
title('Pitch Angle');
legend('Actual', 'Command');
grid on;

subplot(3,2,5);
plot(t, euler_history(3,:), 'b-', 'LineWidth', 1.5);
hold on;
plot(t, euler_cmd_history(3,:), 'r--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Yaw (deg)');
title('Yaw Angle');
legend('Actual', 'Command');
grid on;

subplot(3,2,2);
plot(t, theta_err_history, 'k-', 'LineWidth', 1.5);
hold on;
yline(1, 'g--', '1° threshold');
xline(settling_time_actual, 'm--', 'Settling');
xlabel('Time (s)'); ylabel('Error (deg)');
title('Total Attitude Error');
grid on;

subplot(3,2,4);
plot(t, omega_history(1,:)*180/pi, 'r-', 'LineWidth', 1);
hold on;
plot(t, omega_history(2,:)*180/pi, 'g-', 'LineWidth', 1);
plot(t, omega_history(3,:)*180/pi, 'b-', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Rate (deg/s)');
title('Angular Velocity');
legend('\omega_x', '\omega_y', '\omega_z');
grid on;

subplot(3,2,6);
plot(t, tau_cmd_history(1,:)*1000, 'r-', 'LineWidth', 1);
hold on;
plot(t, tau_cmd_history(2,:)*1000, 'g-', 'LineWidth', 1);
plot(t, tau_cmd_history(3,:)*1000, 'b-', 'LineWidth', 1);
yline(actuators.wheel_max_torque*1000, 'k--');
yline(-actuators.wheel_max_torque*1000, 'k--');
xlabel('Time (s)'); ylabel('Torque (mNm)');
title('Control Torque Commands');
legend('\tau_x', '\tau_y', '\tau_z');
grid on;

saveas(gcf, 'attitude_control_response.png');
fprintf('   Saved: attitude_control_response.png\n');

% Figure 2: Actuator Response
figure('Name', 'Actuator Response', 'Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(t, wheel_cmd_history(1,:)*1000, 'r-', 'LineWidth', 1);
hold on;
plot(t, wheel_cmd_history(2,:)*1000, 'g-', 'LineWidth', 1);
plot(t, wheel_cmd_history(3,:)*1000, 'b-', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Torque (mNm)');
title('Reaction Wheel Torque Commands');
legend('Wheel X', 'Wheel Y', 'Wheel Z');
grid on;

subplot(2,2,2);
plot(t, h_wheels_history(1,:), 'r-', 'LineWidth', 1.5);
hold on;
plot(t, h_wheels_history(2,:), 'g-', 'LineWidth', 1.5);
plot(t, h_wheels_history(3,:), 'b-', 'LineWidth', 1.5);
yline(actuators.wheel_max_momentum, 'k--', 'Max');
yline(-actuators.wheel_max_momentum, 'k--');
xlabel('Time (s)'); ylabel('Momentum (Nms)');
title('Wheel Angular Momentum');
legend('h_x', 'h_y', 'h_z');
grid on;

subplot(2,2,3);
plot(t, mtq_cmd_history(1,:), 'r-', 'LineWidth', 1);
hold on;
plot(t, mtq_cmd_history(2,:), 'g-', 'LineWidth', 1);
plot(t, mtq_cmd_history(3,:), 'b-', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Dipole (Am²)');
title('Magnetorquer Commands');
legend('MTQ X', 'MTQ Y', 'MTQ Z');
grid on;

subplot(2,2,4);
% Sun vector in body frame
sun_body = zeros(3, N);
for k = 1:N
    A = quat_to_dcm(q_history(:,k));
    sun_body(:,k) = A * sun_eci_history(:,k);
end
plot(t, sun_body(1,:), 'r-', 'LineWidth', 1);
hold on;
plot(t, sun_body(2,:), 'g-', 'LineWidth', 1);
plot(t, sun_body(3,:), 'b-', 'LineWidth', 1);
yline(1, 'k--', '+Z target');
xlabel('Time (s)'); ylabel('Component');
title('Sun Vector in Body Frame');
legend('s_x', 's_y', 's_z', 'Location', 'best');
grid on;

saveas(gcf, 'actuator_response.png');
fprintf('   Saved: actuator_response.png\n');

fprintf('\n=== TEST TAMAMLANDI ===\n');

%% ==================== LOCAL HELPER FUNCTIONS ====================

function q = euler_to_quat(e)
    cr=cos(e(1)/2); sr=sin(e(1)/2);
    cp=cos(e(2)/2); sp=sin(e(2)/2);
    cy=cos(e(3)/2); sy=sin(e(3)/2);
    q = [sr*cp*cy-cr*sp*sy; cr*sp*cy+sr*cp*sy; cr*cp*sy-sr*sp*cy; cr*cp*cy+sr*sp*sy];
    q = q/norm(q);
end

function e = quat_to_euler(q)
    q1=q(1); q2=q(2); q3=q(3); q4=q(4);
    e = [atan2(2*(q4*q1+q2*q3), 1-2*(q1^2+q2^2));
         asin(max(-1,min(1,2*(q4*q2-q3*q1))));
         atan2(2*(q4*q3+q1*q2), 1-2*(q2^2+q3^2))];
end

function A = quat_to_dcm(q)
    q1=q(1); q2=q(2); q3=q(3); q4=q(4);
    A = [q4^2+q1^2-q2^2-q3^2, 2*(q1*q2+q3*q4), 2*(q1*q3-q2*q4);
         2*(q1*q2-q3*q4), q4^2-q1^2+q2^2-q3^2, 2*(q2*q3+q1*q4);
         2*(q1*q3+q2*q4), 2*(q2*q3-q1*q4), q4^2-q1^2-q2^2+q3^2];
end

function Omega = quat_omega_matrix(w)
    Omega = [0, w(3), -w(2), w(1);
             -w(3), 0, w(1), w(2);
             w(2), -w(1), 0, w(3);
             -w(1), -w(2), -w(3), 0];
end

function s = compute_sun_position(utc)
    JD = juliandate(utc);
    T = (JD - 2451545) / 36525;
    L = mod(280.46 + 36000.77*T, 360);
    M = mod(357.53 + 35999.05*T, 360) * pi/180;
    lam = (L + 1.915*sin(M) + 0.02*sin(2*M)) * pi/180;
    eps = (23.439 - 0.013*T) * pi/180;
    s = [cos(lam); cos(eps)*sin(lam); sin(eps)*sin(lam)];
    s = s / norm(s);
end

function B = compute_mag_field(r)
    B0 = 3.12e-5;
    Re = 6378137;
    r_mag = norm(r);
    r_hat = r / r_mag;
    m_hat = [0; 0; 1];
    B = B0 * (Re/r_mag)^3 * (3*dot(m_hat, r_hat)*r_hat - m_hat);
end

function [q_cmd, info] = sun_pointing_guidance(sun_eci, ~, ~, ~, params)
    % Simplified sun pointing for test
    body_sun = params.body_axis_to_sun / norm(params.body_axis_to_sun);
    sun_eci = sun_eci / norm(sun_eci);
    
    % Build frame
    z_eci = sun_eci;
    ref_secondary = [0; 0; 1];
    if abs(dot(z_eci, ref_secondary)) > 0.99
        ref_secondary = [0; 1; 0];
    end
    y_eci = cross(z_eci, ref_secondary);
    y_eci = y_eci / norm(y_eci);
    x_eci = cross(y_eci, z_eci);
    
    z_body = body_sun;
    if abs(z_body(3)) < 0.99
        y_body = cross(z_body, [0;0;1]);
    else
        y_body = cross(z_body, [0;1;0]);
    end
    y_body = y_body / norm(y_body);
    x_body = cross(y_body, z_body);
    
    A_eci = [x_eci, y_eci, z_eci]';
    A_body = [x_body, y_body, z_body];
    A_eci2body = A_body' * A_eci;
    
    % DCM to quaternion (Shepperd)
    tr = trace(A_eci2body);
    S = [1+tr; 1+2*A_eci2body(1,1)-tr; 1+2*A_eci2body(2,2)-tr; 1+2*A_eci2body(3,3)-tr];
    [~,idx] = max(S);
    A = A_eci2body;
    switch idx
        case 1
            q4=0.5*sqrt(S(1)); q1=(A(2,3)-A(3,2))/(4*q4); q2=(A(3,1)-A(1,3))/(4*q4); q3=(A(1,2)-A(2,1))/(4*q4);
        case 2
            q1=0.5*sqrt(S(2)); q4=(A(2,3)-A(3,2))/(4*q1); q2=(A(1,2)+A(2,1))/(4*q1); q3=(A(3,1)+A(1,3))/(4*q1);
        case 3
            q2=0.5*sqrt(S(3)); q4=(A(3,1)-A(1,3))/(4*q2); q1=(A(1,2)+A(2,1))/(4*q2); q3=(A(2,3)+A(3,2))/(4*q2);
        case 4
            q3=0.5*sqrt(S(4)); q4=(A(1,2)-A(2,1))/(4*q3); q1=(A(3,1)+A(1,3))/(4*q3); q2=(A(2,3)+A(3,2))/(4*q3);
    end
    q_cmd = [q1;q2;q3;q4];
    if q_cmd(4)<0, q_cmd=-q_cmd; end
    q_cmd = q_cmd/norm(q_cmd);
    info = [];
end

function [tau_cmd, ctrl_info] = pd_attitude_controller(q_cur, q_des, w_cur, w_des, params)
    q_cur = q_cur/norm(q_cur);
    q_des = q_des/norm(q_des);
    
    % Quaternion error
    q_des_conj = [-q_des(1:3); q_des(4)];
    q_err = quat_mult(q_des_conj, q_cur);
    if q_err(4)<0, q_err=-q_err; end
    
    q_err_v = q_err(1:3);
    w_err = w_cur - w_des;
    
    Kp = params.Kp; Kd = params.Kd;
    if isscalar(Kp), Kp=Kp*eye(3); end
    if isscalar(Kd), Kd=Kd*eye(3); end
    
    tau_p = -Kp * sign(q_err(4)) * q_err_v;
    tau_d = -Kd * w_err;
    
    if params.feedforward && isfield(params,'h_wheels')
        L = params.I*w_cur + params.h_wheels;
        tau_ff = cross(w_cur, L);
    else
        tau_ff = [0;0;0];
    end
    
    tau_cmd = tau_p + tau_d + tau_ff;
    
    % Saturate
    tau_max = params.tau_max;
    for i=1:3
        if abs(tau_cmd(i))>tau_max
            tau_cmd(i)=sign(tau_cmd(i))*tau_max;
        end
    end
    
    ctrl_info.theta_error_deg = 2*acos(min(1,abs(q_err(4))))*180/pi;
end

function r = quat_mult(p, q)
    r = [p(4)*q(1)+p(1)*q(4)+p(2)*q(3)-p(3)*q(2);
         p(4)*q(2)-p(1)*q(3)+p(2)*q(4)+p(3)*q(1);
         p(4)*q(3)+p(1)*q(2)-p(2)*q(1)+p(3)*q(4);
         p(4)*q(4)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)];
end

function [wheel_cmd, mtq_cmd, info] = actuator_allocation(tau_cmd, h_wheels, B_body, params)
    wheel_cmd = -tau_cmd;  % Direct allocation
    tau_max = params.wheel_max_torque;
    for i=1:3
        if abs(wheel_cmd(i))>tau_max
            wheel_cmd(i)=sign(wheel_cmd(i))*tau_max;
        end
    end
    
    % Momentum dumping
    mtq_cmd = [0;0;0];
    h_thresh = params.momentum_threshold;
    if any(abs(h_wheels) > h_thresh) && norm(B_body) > 1e-9
        tau_dump = -params.momentum_dump_gain * h_wheels;
        B_norm = norm(B_body);
        tau_perp = tau_dump - dot(tau_dump,B_body)/B_norm^2*B_body;
        mtq_cmd = cross(B_body, tau_perp) / B_norm^2;
        m_max = params.mtq_max_dipole;
        for i=1:3
            if abs(mtq_cmd(i))>m_max
                mtq_cmd(i)=sign(mtq_cmd(i))*m_max;
            end
        end
    end
    info = [];
end