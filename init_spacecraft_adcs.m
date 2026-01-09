%% ========================================================================
%  SPACECRAFT ADCS - INITIALIZATION SCRIPT
%  ========================================================================
%  Bu script, Simulink modeli çalışmadan önce tüm parametreleri yükler.
%  
%  Kullanım:
%    1. Bu dosyayı MATLAB path'ine ekleyin
%    2. Simulink Model Properties > Callbacks > InitFcn: init_spacecraft_adcs
%    3. Veya manuel olarak çalıştırın: >> init_spacecraft_adcs
%
%  Author: Mehmet Yasin Meriç
%  Date: 2026
%% ========================================================================

function init_spacecraft_adcs()

disp('═══════════════════════════════════════════════════════════════════');
disp('   SPACECRAFT ADCS INITIALIZATION');
disp('   60 kg LEO Satellite - Sun Pointing Control');
disp('═══════════════════════════════════════════════════════════════════');

%% ==================== SIMULATION PARAMETERS ====================
sim_params.dt = 0.1;              % Sample time (s) - 10 Hz
sim_params.t_final = 600;         % Simulation duration (s)
sim_params.solver = 'ode4';       % Fixed-step Runge-Kutta

assignin('base', 'sim_params', sim_params);

%% ==================== SPACECRAFT PARAMETERS ====================
spacecraft.mass = 60;                              % kg
spacecraft.I = diag([5.4, 5.4, 4.32]);            % kg·m² (principal)
spacecraft.I_inv = inv(spacecraft.I);
spacecraft.cg = [0; 0; 0];                        % Center of gravity (body)

assignin('base', 'spacecraft', spacecraft);

%% ==================== ORBIT PARAMETERS ====================
orbit.altitude = 1000e3;                           % m
orbit.R_earth = 6378137;                           % m
orbit.mu = 3.986004418e14;                         % m³/s²
orbit.a = orbit.R_earth + orbit.altitude;          % Semi-major axis (m)
orbit.e = 0;                                       % Eccentricity (circular)
orbit.i = 60 * pi/180;                             % Inclination (rad)
orbit.RAAN = 0;                                    % Right ascension (rad)
orbit.omega = 0;                                   % Argument of perigee (rad)
orbit.M0 = 0;                                      % Initial mean anomaly (rad)
orbit.n = sqrt(orbit.mu / orbit.a^3);              % Mean motion (rad/s)
orbit.period = 2*pi / orbit.n;                     % Orbital period (s)

assignin('base', 'orbit', orbit);

%% ==================== INITIAL CONDITIONS ====================
% Attitude (quaternion: ECI to Body) - scalar-last convention
init.euler_deg = [15; -10; 20];                    % Initial Euler angles (deg)
init.euler_rad = init.euler_deg * pi/180;
init.q = euler2quat(init.euler_rad);               % [q1;q2;q3;q4]

% Angular velocity (body frame, rad/s)
init.omega = [0.2; -0.1; 0.15] * pi/180;

% Reaction wheel momentum (Nms)
init.h_wheels = [0; 0; 0];

% Fuel slosh initial state
init.psi = 0;                                      % Slosh angle (rad)
init.psi_dot = 0;                                  % Slosh rate (rad/s)

% MEKF initial state
init.q_est = init.q;                               % Initial estimate = true
init.gyro_bias = [0; 0; 0];                        % Bias estimate
init.P = blkdiag(eye(3)*(10*pi/180)^2, eye(3)*(1*pi/180/3600)^2);

assignin('base', 'init', init);

%% ==================== SENSOR PARAMETERS ====================
% IMU (STIM377H specifications)
sensors.imu.gyro_arw = 0.15 * pi/180 / sqrt(3600);     % rad/s/√Hz
sensors.imu.gyro_bias_instability = 0.5 * pi/180/3600; % rad/s
sensors.imu.gyro_bias_true = [0.4; -0.3; 0.5] * pi/180/3600; % True bias
sensors.imu.accel_noise = 100e-6 * 9.81;               % m/s²
sensors.imu.sample_rate = 100;                          % Hz

% Magnetometer (MAG-3 specifications)
sensors.mag.noise_density = 7e-12;                      % T/√Hz
sensors.mag.noise_std = 10e-9;                          % T (typical)
sensors.mag.sample_rate = 100;                          % Hz

% GNSS (OEM719 specifications)
sensors.gnss.pos_accuracy = 1.5;                        % m (1σ)
sensors.gnss.vel_accuracy = 0.03;                       % m/s (1σ)
sensors.gnss.sample_rate = 1;                           % Hz

assignin('base', 'sensors', sensors);

%% ==================== ACTUATOR PARAMETERS ====================
% Reaction Wheels (3x orthogonal configuration)
actuators.rw.max_torque = 0.01;                % Nm (10 mNm)
actuators.rw.max_momentum = 0.4;               % Nms
actuators.rw.max_speed = 6000 * 2*pi/60;       % rad/s
actuators.rw.inertia = 1e-4;                   % kg·m² (per wheel)
actuators.rw.friction = 1e-6;                  % Nms (viscous)
actuators.rw.num_wheels = 3;
actuators.rw.axes = eye(3);                    % Wheel axes in body frame

% Magnetorquers (3x orthogonal configuration)
actuators.mtq.max_dipole = 5;                  % Am²
actuators.mtq.residual_dipole = [0.1; 0.1; 0.1]; % Am² (residual)
actuators.mtq.num_mtq = 3;
actuators.mtq.axes = eye(3);                   % MTQ axes in body frame

assignin('base', 'actuators', actuators);

%% ==================== CONTROLLER PARAMETERS ====================
% PD Controller gains (torque-limited design)
I_avg = mean(diag(spacecraft.I));
q_err_max = 0.4;  % Max quaternion vector error

controller.Kp = actuators.rw.max_torque / q_err_max * 0.9;
controller.Kd = 2 * 0.707 * sqrt(controller.Kp * I_avg);
controller.Ki = 0;  % No integral term (PD only)
controller.tau_max = actuators.rw.max_torque;
controller.feedforward = 1;  % Enable gyroscopic feedforward

% Actuator allocation
controller.momentum_dump_gain = 0.05;
controller.momentum_threshold = 0.7 * actuators.rw.max_momentum;
controller.allocation_mode = 1;  % 1=hybrid, 2=wheel_only, 3=mtq_only

assignin('base', 'controller', controller);

%% ==================== MEKF PARAMETERS ====================
mekf_params = struct();
mekf_params.gyro_noise = 0.15 * (pi/180) / sqrt(3600);
mekf_params.gyro_bias_noise = 0.5 * (pi/180) / 3600;
mekf_params.accel_noise = 0.03;
mekf_params.mag_noise = 0.05;
mekf_params.init_attitude_var = (10 * pi/180)^2;
mekf_params.init_bias_var = (1 * pi/180 / 3600)^2;
mekf_params.mekf_state = [];

% mekf.gyro_noise = sensors.imu.gyro_arw;
% mekf.gyro_bias_noise = sensors.imu.gyro_bias_instability;
% mekf.accel_noise = 0.02;  % Normalized measurement noise
% mekf.mag_noise = 0.03;    % Normalized measurement noise
% mekf.dt = sim_params.dt;

assignin('base', 'mekf_params', mekf_params);

%% ==================== ENVIRONMENT PARAMETERS ====================
env.B0 = 3.12e-5;              % Earth magnetic field constant (T)
env.Re = 6378137;              % Earth radius (m)
env.mu = orbit.mu;             % Gravitational parameter
env.P_sun = 4.56e-6;           % Solar pressure at 1 AU (N/m²)

assignin('base', 'env', env);

%% ==================== DISTURBANCE PARAMETERS ====================
disturbance.gravity_gradient = 1;      % Enable (1) / Disable (0)
disturbance.magnetic = 1;
disturbance.solar_pressure = 1;
disturbance.aerodynamic = 0;           % Negligible at 1000 km

assignin('base', 'disturbance', disturbance);

%% ==================== FUEL SLOSH PARAMETERS ====================
slosh.enabled = 1;                     % Enable (1) / Disable (0)
slosh.m_slosh = 2.6;                   % Sloshing mass (kg)
slosh.m_fixed = 1.4;                   % Fixed fuel mass (kg)
slosh.L = 0.084;                       % Pendulum length (m)
slosh.zeta = 0.01;                     % Damping ratio
slosh.tank_position = [0; 0; 0.15];    % Tank position from CoM (m)
slosh.g_eff = 0.05 * 9.81;             % Effective gravity (m/s²)

assignin('base', 'slosh', slosh);

%% ==================== REFERENCE TIME ====================
ref_time.year = 2026;
ref_time.month = 6;
ref_time.day = 21;
ref_time.hour = 12;
ref_time.min = 0;
ref_time.sec = 0;
ref_time.jd0 = juliandate(datetime(ref_time.year, ref_time.month, ref_time.day, ...
                                   ref_time.hour, ref_time.min, ref_time.sec));

assignin('base', 'ref_time', ref_time);

%% ==================== POINTING MODE ====================
guidance.mode = 1;                     % 1=sun, 2=nadir, 3=inertial
guidance.body_axis_to_sun = [0; 0; 1]; % +Z axis to sun

assignin('base', 'guidance', guidance);

%% ==================== DISPLAY SUMMARY ====================
disp('───────────────────────────────────────────────────────────────────');
disp('  CONFIGURATION SUMMARY');
disp('───────────────────────────────────────────────────────────────────');
fprintf('  Spacecraft mass:       %.1f kg\n', spacecraft.mass);
fprintf('  Inertia (diag):        [%.2f, %.2f, %.2f] kg·m²\n', ...
        spacecraft.I(1,1), spacecraft.I(2,2), spacecraft.I(3,3));
fprintf('  Orbit altitude:        %.0f km\n', orbit.altitude/1000);
fprintf('  Orbit period:          %.1f min\n', orbit.period/60);
fprintf('  Orbit inclination:     %.1f deg\n', orbit.i*180/pi);
disp('───────────────────────────────────────────────────────────────────');
fprintf('  Initial attitude:      [%.1f, %.1f, %.1f] deg\n', init.euler_deg);
fprintf('  Initial rate:          [%.2f, %.2f, %.2f] deg/s\n', init.omega*180/pi);
disp('───────────────────────────────────────────────────────────────────');
fprintf('  RW max torque:         %.1f mNm\n', actuators.rw.max_torque*1000);
fprintf('  RW max momentum:       %.2f Nms\n', actuators.rw.max_momentum);
fprintf('  MTQ max dipole:        %.1f Am²\n', actuators.mtq.max_dipole);
disp('───────────────────────────────────────────────────────────────────');
fprintf('  Controller Kp:         %.4f Nm/rad\n', controller.Kp);
fprintf('  Controller Kd:         %.4f Nm/(rad/s)\n', controller.Kd);
disp('───────────────────────────────────────────────────────────────────');
fprintf('  Simulation duration:   %.0f s (%.1f min)\n', sim_params.t_final, sim_params.t_final/60);
fprintf('  Sample time:           %.3f s (%.0f Hz)\n', sim_params.dt, 1/sim_params.dt);
disp('═══════════════════════════════════════════════════════════════════');
disp('  INITIALIZATION COMPLETE');
disp('═══════════════════════════════════════════════════════════════════');

end

%% ==================== HELPER FUNCTIONS ====================

function q = euler2quat(e)
% EULER2QUAT - Convert Euler angles (3-2-1) to quaternion
% Input:  e = [roll; pitch; yaw] in radians
% Output: q = [q1; q2; q3; q4] scalar-last convention

    cr = cos(e(1)/2); sr = sin(e(1)/2);
    cp = cos(e(2)/2); sp = sin(e(2)/2);
    cy = cos(e(3)/2); sy = sin(e(3)/2);
    
    q = [sr*cp*cy - cr*sp*sy;
         cr*sp*cy + sr*cp*sy;
         cr*cp*sy - sr*sp*cy;
         cr*cp*cy + sr*sp*sy];
    
    q = q / norm(q);
    
    % Ensure scalar part is positive
    if q(4) < 0
        q = -q;
    end
end
