% Fuel Slosh Model - Pendulum Analog
% Case Study (f): Uzay aracındaki tankın kütlesini ve tankın içindeki 
% sıvıyı belirleyiniz ve bu sıvının sistem üzerindeki etkisini modelleyiniz.
%
% Bu model, yakıt çalkalanmasını (fuel slosh) pendulum analog ile simüle eder.
% Spherical tank içindeki sıvının first-mode slosh davranışı modellenir.
%
% Referanslar:
%   - Dodge, F.T., "The New Dynamic Behavior of Liquids in Moving Containers"
%   - NASA SP-106, "Propellant Slosh Loads"
%   - Abramson, H.N., "The Dynamic Behavior of Liquids in Moving Containers"
%   - MathWorks, "Modeling Spacecraft Fuel Slosh at Embry-Riddle"
%
% Mekanik Analog Modeli:
%   - Sloshing sıvı kütlesi → Pendulum bob mass
%   - Slosh frekansı → Pendulum doğal frekansı
%   - Viskoz sönümleme → Dashpot damping
%
% Denklemler:
%   Pendulum açısal hareketi (küçük açı):
%   (m_slosh * L^2) * ψ̈ + c * ψ̇ + m_slosh * g_eff * L * ψ = -m_slosh * L * θ̈
%   
%   Spacecraft'a etki eden tork:
%   τ_slosh = -m_slosh * L^2 * (ψ̈ + θ̈) - m_slosh * g_eff * L * sin(ψ)

%% ==================== TANK VE SIVI TANIMLARI ====================
%
% TASARIM KARARLARI:
%
% Tank Tipi: Spherical (küresel)
%   - Avantaj: Her yönde simetrik, slosh analizi daha basit
%   - Dezavantaj: Volume utilization düşük
%
% Sıvı: Hydrazine (N2H4) - Monopropellant
%   - Yoğunluk: 1004 kg/m³ @ 25°C
%   - Kinematik viskozite: 0.97 cSt (9.7×10⁻⁷ m²/s)
%   - Yaygın kullanım: Attitude control thrusters
%   - Alternatifler: MMH, NTO, Xenon (ion propulsion)
%
% Tank Doluluk: %50 (case study'de belirtildiği gibi)
%

function [slosh_params, tank_params] = fuel_slosh_tank_design(spacecraft_mass)
% FUEL_SLOSH_TANK_DESIGN - Tank ve sıvı parametrelerini tasarlar
%
% Girdiler:
%   spacecraft_mass - toplam uzay aracı kütlesi (kg) [default: 60]
%
% Çıktılar:
%   slosh_params - slosh model parametreleri (struct)
%   tank_params  - tank fiziksel parametreleri (struct)

if nargin < 1
    spacecraft_mass = 60;  % kg (case study)
end

fprintf('=================================================================\n');
fprintf('   FUEL SLOSH TANK TASARIMI\n');
fprintf('   Case Study (f): Tank kütlesi ve sıvı etkileri\n');
fprintf('=================================================================\n\n');

%% ==================== PROPELLANT SEÇİMİ ====================
fprintf('1. Propellant Seçimi: HYDRAZINE (N2H4)\n');

% Hydrazine özellikleri
propellant.name = 'Hydrazine (N2H4)';
propellant.density = 1004;          % kg/m³ @ 25°C
propellant.kinematic_viscosity = 9.7e-7;  % m²/s (0.97 cSt)
propellant.surface_tension = 0.0663;  % N/m
propellant.specific_impulse = 220;    % s (monopropellant)

fprintf('   - Yoğunluk: %.0f kg/m³\n', propellant.density);
fprintf('   - Kinematik viskozite: %.2f cSt\n', propellant.kinematic_viscosity*1e6);
fprintf('   - Yüzey gerilimi: %.4f N/m\n', propellant.surface_tension);
fprintf('   - Özgül itki (Isp): %.0f s\n', propellant.specific_impulse);

%% ==================== TANK BOYUTLANDIRMA ====================
fprintf('\n2. Tank Boyutlandırma\n');

% Tipik small satellite propellant budget: %10-20 of dry mass
% 60 kg spacecraft için ~6-10 kg propellant
propellant_mass_total = 8;  % kg (conservative estimate)
fill_ratio = 0.50;          % %50 dolu (case study requirement)

% Current propellant mass at 50% fill
propellant_mass = propellant_mass_total * fill_ratio;  % kg

% Tank volume calculation
propellant_volume = propellant_mass / propellant.density;  % m³

% Spherical tank for simplicity
% Total tank volume (account for ullage at 50% fill)
tank_volume_total = propellant_volume / fill_ratio;  % m³

% Tank radius
tank_radius = (3 * tank_volume_total / (4 * pi))^(1/3);  % m

% Tank mass (Titanium, typical wall thickness ~1mm for this size)
tank_wall_thickness = 0.001;  % m (1 mm)
titanium_density = 4500;      % kg/m³
tank_surface_area = 4 * pi * tank_radius^2;
tank_mass = tank_surface_area * tank_wall_thickness * titanium_density;

fprintf('   - Toplam propellant kapasitesi: %.1f kg\n', propellant_mass_total);
fprintf('   - Mevcut propellant (50%%): %.1f kg\n', propellant_mass);
fprintf('   - Propellant hacmi: %.4f m³ (%.1f L)\n', propellant_volume, propellant_volume*1000);
fprintf('   - Tank iç yarıçapı: %.3f m (%.1f cm)\n', tank_radius, tank_radius*100);
fprintf('   - Tank kütlesi (Ti): %.2f kg\n', tank_mass);

% Store tank parameters
tank_params.propellant = propellant;
tank_params.radius = tank_radius;
tank_params.volume_total = tank_volume_total;
tank_params.volume_liquid = propellant_volume;
tank_params.fill_ratio = fill_ratio;
tank_params.tank_mass = tank_mass;
tank_params.propellant_mass = propellant_mass;
tank_params.wall_thickness = tank_wall_thickness;

%% ==================== SLOSH PARAMETRELERİ (PENDULUM ANALOG) ====================
fprintf('\n3. Slosh Parametreleri (Pendulum Analog)\n');

% Dodge formülasyonu (spherical tank, lateral slosh, first mode)
% Reference: NASA SP-106, Dodge "Dynamic Behavior of Liquids"

% Fill ratio dependent parameters for spherical tank
% h = liquid height measured from tank bottom
h_liquid = 2 * tank_radius * fill_ratio;  % liquid depth

% Non-dimensional fill depth
h_over_d = h_liquid / (2 * tank_radius);  % = fill_ratio for sphere

% First slosh mode natural frequency (Dodge, spherical tank)
% ω_n = sqrt(g_eff / L_eff)
% For 50% fill spherical tank: ω_n ≈ sqrt(1.19 * g / R)
% Bu coefficient fill ratio'ya bağlı değişir

% Empirical coefficient for spherical tank (from Dodge tables)
% At 50% fill: k_freq ≈ 1.19
k_freq = 1.19;

% Effective pendulum length
% L_eff = R / k_freq (yaklaşık)
L_pendulum = tank_radius / k_freq;

% Slosh mass (participating mass in first mode)
% For spherical tank at 50% fill: m_slosh ≈ 0.65 * m_liquid
k_mass = 0.65;  % mass participation factor
m_slosh = k_mass * propellant_mass;

% Fixed (non-sloshing) mass
m_fixed = propellant_mass - m_slosh;

% Pendulum hinge location (below liquid surface)
% For spherical tank at 50% fill: h_hinge ≈ 0.75 * R from tank center
k_hinge = 0.75;
h_hinge = k_hinge * tank_radius;

% Damping ratio estimation (depends on viscosity and tank geometry)
% For smooth tank walls: ζ ≈ 0.5-2% (very low)
% With baffles: ζ ≈ 5-15%
% We assume smooth tank (no baffles)
damping_ratio = 0.01;  % 1% (typical for smooth spherical tank)

% Natural frequency at 1g (ground test reference)
g_earth = 9.81;
omega_n_1g = sqrt(k_freq * g_earth / tank_radius);
freq_n_1g = omega_n_1g / (2*pi);

% Slosh in microgravity (different behavior)
% For typical LEO maneuver accelerations: 0.01-0.1 g
g_eff_typical = 0.05 * g_earth;  % 50 milli-g
omega_n_orbit = sqrt(k_freq * g_eff_typical / tank_radius);
freq_n_orbit = omega_n_orbit / (2*pi);

fprintf('   - Slosh kütlesi (m_slosh): %.2f kg (%.0f%% of liquid)\n', ...
        m_slosh, k_mass*100);
fprintf('   - Sabit kütle (m_fixed): %.2f kg\n', m_fixed);
fprintf('   - Pendulum uzunluğu (L): %.4f m (%.1f cm)\n', ...
        L_pendulum, L_pendulum*100);
fprintf('   - Hinge pozisyonu: %.4f m (tank merkezinden)\n', h_hinge);
fprintf('   - Sönüm oranı (ζ): %.1f%%\n', damping_ratio*100);
fprintf('   - Doğal frekans @1g: %.3f Hz (%.2f rad/s)\n', freq_n_1g, omega_n_1g);
fprintf('   - Doğal frekans @0.05g: %.3f Hz (%.2f rad/s)\n', freq_n_orbit, omega_n_orbit);

% Store slosh parameters
slosh_params.m_slosh = m_slosh;           % Sloshing mass (kg)
slosh_params.m_fixed = m_fixed;           % Fixed mass (kg)
slosh_params.L = L_pendulum;              % Pendulum length (m)
slosh_params.h_hinge = h_hinge;           % Hinge location from center (m)
slosh_params.zeta = damping_ratio;        % Damping ratio
slosh_params.k_freq = k_freq;             % Frequency coefficient
slosh_params.k_mass = k_mass;             % Mass participation factor
slosh_params.omega_n_1g = omega_n_1g;     % Natural frequency at 1g
slosh_params.omega_n = omega_n_orbit;     % Natural frequency at typical orbit accel

% Damping coefficient
slosh_params.c = 2 * damping_ratio * omega_n_orbit * m_slosh * L_pendulum^2;

%% ==================== KÜTLE ÖZELLİKLERİ ETKİSİ ====================
fprintf('\n4. Kütle Özellikleri Üzerindeki Etki\n');

% Tank pozisyonu (spacecraft body frame)
% Assume tank is on +Z axis (typical for thrusters)
tank_position = [0; 0; 0.15];  % m (15 cm from CoM)

% Spacecraft inertia without fuel (rough estimate)
% Typical small sat: I ≈ m * r² where r ≈ 0.3m (cube dimension)
I_sc_dry = spacecraft_mass * (0.3)^2 * diag([1, 1, 0.8]);  % kg·m²

% Inertia contribution from tank (fixed mass at tank center)
r_tank = norm(tank_position);
I_tank_fixed = m_fixed * (eye(3) * r_tank^2 - tank_position * tank_position');

% Inertia contribution from slosh mass (at pendulum end)
% This varies with slosh angle ψ
% For small ψ: I_slosh ≈ m_slosh * (r_tank + L)²
I_slosh_nominal = m_slosh * (r_tank + L_pendulum)^2 * eye(3);

% Total inertia
I_total = I_sc_dry + I_tank_fixed + I_slosh_nominal;

fprintf('   - Tank pozisyonu: [%.2f, %.2f, %.2f] m\n', tank_position);
fprintf('   - S/C kuru ataletleri: [%.3f, %.3f, %.3f] kg·m²\n', ...
        I_sc_dry(1,1), I_sc_dry(2,2), I_sc_dry(3,3));
fprintf('   - Tank+Fuel ataleti: [%.3f, %.3f, %.3f] kg·m²\n', ...
        I_tank_fixed(1,1)+I_slosh_nominal(1,1), ...
        I_tank_fixed(2,2)+I_slosh_nominal(2,2), ...
        I_tank_fixed(3,3)+I_slosh_nominal(3,3));
fprintf('   - Toplam atalet: [%.3f, %.3f, %.3f] kg·m²\n', ...
        I_total(1,1), I_total(2,2), I_total(3,3));

% Store mass properties
slosh_params.tank_position = tank_position;
slosh_params.I_sc_dry = I_sc_dry;
slosh_params.I_total_nominal = I_total;

%% ==================== ÖZET ====================
fprintf('\n=================================================================\n');
fprintf('   TASARIM ÖZETİ\n');
fprintf('=================================================================\n');
fprintf('   %-25s: %s\n', 'Propellant', propellant.name);
fprintf('   %-25s: %.1f L\n', 'Tank hacmi', tank_volume_total*1000);
fprintf('   %-25s: %.1f cm\n', 'Tank çapı', 2*tank_radius*100);
fprintf('   %-25s: %.2f kg\n', 'Tank kütlesi', tank_mass);
fprintf('   %-25s: %.2f kg\n', 'Propellant kütlesi (50%%)', propellant_mass);
fprintf('   %-25s: %.2f kg\n', 'Slosh kütlesi', m_slosh);
fprintf('   %-25s: %.3f Hz\n', 'Slosh frekansı (@0.05g)', freq_n_orbit);
fprintf('   %-25s: %.1f%%\n', 'Sönüm oranı', damping_ratio*100);
fprintf('=================================================================\n');

end


%% ==================== SLOSH DYNAMICS FUNCTION ====================

function [psi_dot, psi_ddot, tau_slosh] = fuel_slosh_dynamics(psi, psi_dot_in, theta_ddot, g_eff, slosh_params)
% FUEL_SLOSH_DYNAMICS - Pendulum slosh dinamiklerini hesaplar
%
% Girdiler:
%   psi        - slosh açısı (rad)
%   psi_dot_in - slosh açısal hızı (rad/s)
%   theta_ddot - spacecraft açısal ivmesi (rad/s²)
%   g_eff      - efektif yerçekimi ivmesi (m/s²)
%   slosh_params - slosh parametreleri struct
%
% Çıktılar:
%   psi_dot    - slosh açısal hızı (rad/s)
%   psi_ddot   - slosh açısal ivmesi (rad/s²)
%   tau_slosh  - spacecraft'a etki eden slosh torku (Nm)
%
% Pendulum denklemi (küçük açı yaklaşımı):
%   I_p * ψ̈ + c * ψ̇ + m*g*L*ψ = -m*L²*θ̈
%
% Nonlinear form:
%   I_p * ψ̈ + c * ψ̇ + m*g*L*sin(ψ) = -m*L²*θ̈*cos(ψ)

% Extract parameters
m = slosh_params.m_slosh;
L = slosh_params.L;
c = slosh_params.c;

% Moment of inertia of pendulum about hinge
I_p = m * L^2;

% Update natural frequency for current g_eff
omega_n_sq = slosh_params.k_freq * g_eff / (slosh_params.L / slosh_params.k_freq);

% Pendulum dynamics (nonlinear)
% I_p * psi_ddot = -c * psi_dot - m*g_eff*L*sin(psi) - m*L^2*theta_ddot*cos(psi)

psi_ddot = (-c * psi_dot_in - m * g_eff * L * sin(psi) - m * L^2 * theta_ddot * cos(psi)) / I_p;

% Output psi_dot unchanged
psi_dot = psi_dot_in;

% Slosh torque on spacecraft (reaction to pendulum motion)
% τ = -m*L² * (ψ̈ + θ̈) - m*g_eff*L*sin(ψ)
tau_slosh = -m * L^2 * (psi_ddot + theta_ddot) - m * g_eff * L * sin(psi);

end


%% ==================== COUPLED DYNAMICS FUNCTION ====================

function [x_dot] = spacecraft_slosh_dynamics(t, x, params, tau_control)
% SPACECRAFT_SLOSH_DYNAMICS - Coupled spacecraft + slosh dynamics
%
% State vector x = [θ; θ_dot; ψ; ψ_dot] (2D planar case)
%   θ: spacecraft attitude angle
%   ψ: slosh pendulum angle
%
% Coupled equations of motion:
%   (I_sc + m_f*b²) * θ̈ - m_s*a*b*cos(ψ) * ψ̈ + m_s*a*b*sin(ψ)*(θ̇+ψ̇)² = τ_control + τ_dist
%   (I_f + m_s*a²) * ψ̈ - m_s*a*b*cos(ψ) * θ̈ + m_s*a*g*sin(ψ) + c*ψ̇ = 0

% Extract state
theta = x(1);
theta_dot = x(2);
psi = x(3);
psi_dot = x(4);

% Parameters
I_sc = params.I_sc;           % Spacecraft inertia (without fuel)
I_f = params.I_f;             % Fuel moment of inertia
m_s = params.slosh.m_slosh;   % Slosh mass
m_f = params.slosh.m_fixed;   % Fixed fuel mass
a = params.slosh.L;           % Pendulum length
b = norm(params.slosh.tank_position);  % Distance to tank
c = params.slosh.c;           % Damping coefficient
g = params.g_eff;             % Effective gravity

% Mass matrix elements
M11 = I_sc + m_f * b^2 + m_s * (b^2 - a*b*cos(psi));
M12 = -m_s * a * b * cos(psi);
M21 = M12;
M22 = I_f + m_s * a^2;

% Coriolis/centrifugal and gravity terms
C1 = m_s * a * b * (theta_dot + psi_dot)^2 * sin(psi);
C2 = m_s * a * g * sin(psi) + c * psi_dot;

% Right-hand side
if nargin < 4 || isempty(tau_control)
    tau_control = 0;
end
RHS = [tau_control + C1; -C2];

% Solve for accelerations
M = [M11, M12; M21, M22];
acc = M \ RHS;

% State derivative
x_dot = [theta_dot; acc(1); psi_dot; acc(2)];

end
