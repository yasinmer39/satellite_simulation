% =========================================================================
% Yönelim Dinamiği Test Scripti
% 60 kg LEO Uydu
% =========================================================================
% Bu script, tüm yönelim dinamiği fonksiyonlarını test eder:
%   1. Euler Dinamik Denklemleri
%   2. Quaternion Kinematiği
%   3. DCM Kinematiği
%   4. Euler Açıları Kinematiği
%   5. Koordinat Dönüşümleri
%   6. Torque-Free Motion
%   7. Attitude Propagation
% =========================================================================

clear; clc; close all;

fprintf('========================================\n');
fprintf('   YÖNELİM DİNAMİĞİ TESTİ\n');
fprintf('   Wertz Chapter 16 Denklemleri\n');
fprintf('========================================\n\n');

%% Uydu Parametreleri
% Asal atalet momentleri (kg·m²)
I1 = 5.0;   % x ekseni
I2 = 4.0;   % y ekseni
I3 = 3.0;   % z ekseni (en küçük - tipik spin ekseni)
I = [I1; I2; I3];
I_matrix = diag(I);

fprintf('Uydu Atalet Momentleri:\n');
fprintf('  I1 = %.1f kg·m², I2 = %.1f kg·m², I3 = %.1f kg·m²\n\n', I1, I2, I3);

%% TEST 1: Euler Dinamik Denklemleri
fprintf('TEST 1: EULER DİNAMİK DENKLEMLERİ\n');
fprintf('----------------------------------\n');

omega = [0.1; 0.05; 0.5];   % rad/s
N = [0.001; 0; 0];          % N·m

omega_dot = euler_equations(omega, I, N);

fprintf('Açısal hız ω = [%.3f, %.3f, %.3f] rad/s\n', omega(1), omega(2), omega(3));
fprintf('Dış tork N = [%.4f, %.4f, %.4f] N·m\n', N(1), N(2), N(3));
fprintf('Açısal ivme ω̇ = [%.6f, %.6f, %.6f] rad/s²\n\n', ...
    omega_dot(1), omega_dot(2), omega_dot(3));

%% TEST 2: Quaternion Kinematiği
fprintf('TEST 2: QUATERNION KİNEMATİĞİ\n');
fprintf('------------------------------\n');

q = [0; 0; 0; 1];   % identity quaternion
omega_test = [0.1; 0; 0];

q_dot = quaternion_kinematics(q, omega_test);

fprintf('Quaternion q = [%.3f, %.3f, %.3f, %.3f]\n', q(1), q(2), q(3), q(4));
fprintf('Açısal hız ω = [%.3f, %.3f, %.3f] rad/s\n', omega_test(1), omega_test(2), omega_test(3));
fprintf('Quaternion türevi q̇ = [%.6f, %.6f, %.6f, %.6f]\n\n', ...
    q_dot(1), q_dot(2), q_dot(3), q_dot(4));

%% TEST 3: Koordinat Dönüşümleri
fprintf('TEST 3: KOORDİNAT DÖNÜŞÜMLERİ\n');
fprintf('------------------------------\n');

% Test Euler açıları (3-2-1)
euler_test = [30; 20; 10] * pi/180;  % rad

% Euler -> DCM
A_test = euler2dcm_321(euler_test);

% DCM -> Euler
euler_back = dcm2euler_321(A_test);

% Euler -> Quaternion
q_from_euler = euler2quat_321(euler_test);

% Quaternion -> DCM
A_from_quat = quat2dcm_wertz(q_from_euler);

% Quaternion -> Euler
euler_from_quat = quat2euler_321(q_from_euler);

fprintf('Orijinal Euler açıları: [%.1f, %.1f, %.1f] derece\n', ...
    euler_test(1)*180/pi, euler_test(2)*180/pi, euler_test(3)*180/pi);
fprintf('DCM -> Euler: [%.1f, %.1f, %.1f] derece\n', ...
    euler_back(1)*180/pi, euler_back(2)*180/pi, euler_back(3)*180/pi);
fprintf('Quaternion: [%.4f, %.4f, %.4f, %.4f]\n', ...
    q_from_euler(1), q_from_euler(2), q_from_euler(3), q_from_euler(4));
fprintf('Quat -> Euler: [%.1f, %.1f, %.1f] derece\n\n', ...
    euler_from_quat(1)*180/pi, euler_from_quat(2)*180/pi, euler_from_quat(3)*180/pi);

% DCM orthogonalite kontrolü
fprintf('DCM orthogonalite kontrolü:\n');
fprintf('  det(A) = %.6f (olması gereken: 1)\n', det(A_test));
fprintf('  ||A*A'' - I|| = %.2e (olması gereken: ~0)\n\n', norm(A_test*A_test' - eye(3)));

%% TEST 4: Quaternion İşlemleri
fprintf('TEST 4: QUATERNION İŞLEMLERİ\n');
fprintf('-----------------------------\n');

q1 = euler2quat_321([45; 0; 0] * pi/180);
q2 = euler2quat_321([0; 30; 0] * pi/180);

% Çarpım
q_product = quat_multiply(q1, q2);

% Ters
q1_inv = quat_inverse(q1);

% q1 * q1_inv = identity?
q_identity = quat_multiply(q1, q1_inv);

fprintf('q1 = 45° roll, q2 = 30° pitch\n');
fprintf('q1 ⊗ q2 = [%.4f, %.4f, %.4f, %.4f]\n', ...
    q_product(1), q_product(2), q_product(3), q_product(4));
fprintf('q1 ⊗ q1⁻¹ = [%.4f, %.4f, %.4f, %.4f] (identity olmalı)\n\n', ...
    q_identity(1), q_identity(2), q_identity(3), q_identity(4));

%% TEST 5: Torque-Free Motion (Eksenel Simetrik)
fprintf('TEST 5: TORQUE-FREE MOTION\n');
fprintf('---------------------------\n');

% Eksenel simetrik atalet
I_sym = [5.0; 5.0; 3.0];  % I1 = I2

omega0 = [0.1; 0.05; 0.5];  % rad/s
t_sim = linspace(0, 60, 1000);  % 60 saniye

[omega_history, omega_p] = torque_free_motion(t_sim, omega0, I_sym);

fprintf('Başlangıç ω = [%.3f, %.3f, %.3f] rad/s\n', omega0(1), omega0(2), omega0(3));
fprintf('Body nutation rate ωp = %.4f rad/s (%.2f derece/s)\n', omega_p, omega_p*180/pi);
fprintf('Nutation periyodu = %.2f s\n\n', 2*pi/abs(omega_p));

%% TEST 6: Tam Attitude Propagation (ODE)
fprintf('TEST 6: ATTITUDE PROPAGATION (ODE)\n');
fprintf('-----------------------------------\n');

% Başlangıç durumu
q0 = [0; 0; 0; 1];          % identity
omega0_ode = [0.1; 0.05; 0.5];  % rad/s
state0 = [q0; omega0_ode];

% Simülasyon parametreleri
t_span = [0 60];

% ODE çözümü (torksuz)
N_zero = [0; 0; 0];
ode_func = @(t, state) attitude_dynamics_ode(t, state, I, N_zero);

% ode45 ile çöz
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
[t_ode, state_ode] = ode45(ode_func, t_span, state0, options);

fprintf('Simülasyon süresi: %.1f s\n', t_span(2));
fprintf('Adım sayısı: %d\n', length(t_ode));

% Son durumu kontrol et
q_final = state_ode(end, 1:4)';
omega_final = state_ode(end, 5:7)';

fprintf('Son quaternion: [%.4f, %.4f, %.4f, %.4f]\n', ...
    q_final(1), q_final(2), q_final(3), q_final(4));
fprintf('|q| = %.6f (olması gereken: 1)\n', norm(q_final));

% Enerji ve momentum korunumu kontrolü
L0 = angular_momentum(omega0_ode, I);
L_final = angular_momentum(omega_final, I);
E0 = kinetic_energy(omega0_ode, I);
E_final = kinetic_energy(omega_final, I);

fprintf('\nKorunum kontrolü (torksuz hareket):\n');
fprintf('  |L|: %.6f -> %.6f (fark: %.2e)\n', norm(L0), norm(L_final), abs(norm(L0)-norm(L_final)));
fprintf('  Ek:  %.6f -> %.6f (fark: %.2e)\n\n', E0, E_final, abs(E0-E_final));

%% TEST 7: Quaternion Discrete Propagation
fprintf('TEST 7: QUATERNION DISCRETE PROPAGATION\n');
fprintf('----------------------------------------\n');

q_discrete = [0; 0; 0; 1];
omega_const = [0; 0; 0.1];  % Sadece z ekseni etrafında dönüş
dt = 0.1;

fprintf('Sabit ω = [0, 0, 0.1] rad/s, dt = %.2f s\n', dt);
fprintf('Beklenen 10 saniye sonra: θ = %.1f derece\n', 1.0 * 180/pi);

for i = 1:100
    q_discrete = quat_propagate(q_discrete, omega_const, dt);
end

euler_after = quat2euler_321(q_discrete);
fprintf('10 s sonra Euler açıları: [%.2f, %.2f, %.2f] derece\n\n', ...
    euler_after(1)*180/pi, euler_after(2)*180/pi, euler_after(3)*180/pi);

%% Grafikler
fprintf('GRAFİKLER OLUŞTURULUYOR...\n');

figure('Name', 'Yönelim Dinamiği Test Sonuçları', 'Position', [100 100 1200 800]);

% 1. Torque-Free Motion - Açısal Hız
subplot(2,3,1);
plot(t_sim, omega_history(1,:), 'r-', 'LineWidth', 1.5);
hold on;
plot(t_sim, omega_history(2,:), 'g-', 'LineWidth', 1.5);
plot(t_sim, omega_history(3,:), 'b-', 'LineWidth', 1.5);
xlabel('Zaman (s)');
ylabel('Açısal Hız (rad/s)');
title('Torque-Free Motion');
legend('ω_1', 'ω_2', 'ω_3', 'Location', 'best');
grid on;

% 2. Torque-Free Motion - ω1-ω2 düzlemi
subplot(2,3,2);
plot(omega_history(1,:), omega_history(2,:), 'b-', 'LineWidth', 1.5);
xlabel('ω_1 (rad/s)');
ylabel('ω_2 (rad/s)');
title('Body Nutation (ω_1-ω_2 düzlemi)');
axis equal;
grid on;

% 3. ODE Propagation - Quaternion
subplot(2,3,3);
plot(t_ode, state_ode(:,1), 'r-', 'LineWidth', 1);
hold on;
plot(t_ode, state_ode(:,2), 'g-', 'LineWidth', 1);
plot(t_ode, state_ode(:,3), 'b-', 'LineWidth', 1);
plot(t_ode, state_ode(:,4), 'k-', 'LineWidth', 1);
xlabel('Zaman (s)');
ylabel('Quaternion Bileşenleri');
title('Attitude Propagation (ODE)');
legend('q_1', 'q_2', 'q_3', 'q_4', 'Location', 'best');
grid on;

% 4. ODE Propagation - Açısal Hız
subplot(2,3,4);
plot(t_ode, state_ode(:,5), 'r-', 'LineWidth', 1.5);
hold on;
plot(t_ode, state_ode(:,6), 'g-', 'LineWidth', 1.5);
plot(t_ode, state_ode(:,7), 'b-', 'LineWidth', 1.5);
xlabel('Zaman (s)');
ylabel('Açısal Hız (rad/s)');
title('Açısal Hız Tarihi (ODE)');
legend('ω_1', 'ω_2', 'ω_3', 'Location', 'best');
grid on;

% 5. Enerji ve Momentum Korunumu
subplot(2,3,5);
L_history = zeros(length(t_ode), 1);
E_history = zeros(length(t_ode), 1);
for i = 1:length(t_ode)
    omega_i = state_ode(i, 5:7)';
    L_history(i) = norm(angular_momentum(omega_i, I));
    E_history(i) = kinetic_energy(omega_i, I);
end
yyaxis left;
plot(t_ode, L_history, 'b-', 'LineWidth', 1.5);
ylabel('|L| (kg·m²/s)');
yyaxis right;
plot(t_ode, E_history, 'r-', 'LineWidth', 1.5);
ylabel('E_k (J)');
xlabel('Zaman (s)');
title('Korunum Miktarları');
grid on;

% 6. Euler Açıları Zaman Serisi
subplot(2,3,6);
euler_history = zeros(length(t_ode), 3);
for i = 1:length(t_ode)
    q_i = state_ode(i, 1:4)';
    euler_history(i,:) = quat2euler_321(q_i)' * 180/pi;
end
plot(t_ode, euler_history(:,1), 'r-', 'LineWidth', 1.5);
hold on;
plot(t_ode, euler_history(:,2), 'g-', 'LineWidth', 1.5);
plot(t_ode, euler_history(:,3), 'b-', 'LineWidth', 1.5);
xlabel('Zaman (s)');
ylabel('Euler Açıları (derece)');
title('Attitude Euler Açıları (3-2-1)');
legend('φ (roll)', 'θ (pitch)', 'ψ (yaw)', 'Location', 'best');
grid on;

fprintf('\nTüm testler tamamlandı.\n');
fprintf('========================================\n');
