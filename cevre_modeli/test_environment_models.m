% =========================================================================
% Çevre Modeli Test Scripti
% 60 kg LEO Uydu (1000 km yükseklik)
% =========================================================================
% Bu script, tüm çevre modeli fonksiyonlarını test eder:
%   1. Manyetik Alan (Dipol ve IGRF)
%   2. Güneş Pozisyonu
%   3. Güneş Radyasyon Basıncı (SRP)
%   4. Gölge (Eclipse) Hesaplama
%   5. Atmosferik Sürüklenme (Drag)
%   6. Kütleçekimsel Gradyent Torku
% =========================================================================

clear; clc; close all;

%% Sabitler
global mu
mu = 398600;            % Dünya kütleçekimsel parametresi (km³/s²)
R_earth = 6378.137;     % Dünya yarıçapı (km)

%% Uydu Parametreleri
m = 60;                 % Kütle (kg)
A_cross = 0.5;          % Çapraz kesit alanı (m²)
Cd = 2.2;               % Drag katsayısı
Cr = 1.5;               % Reflectivity katsayısı

% Atalet matrisi (kg·m²) - örnek değerler
Ix = 5.0; Iy = 4.0; Iz = 3.0;
I = diag([Ix, Iy, Iz]);

%% Yörünge Parametreleri
h_altitude = 1000;      % km
a = R_earth + h_altitude;
incl = 98 * pi/180;     % rad (güneş-senkron)

% Başlangıç pozisyonu ve hızı (dairesel yörünge, equator üzerinde)
v_orbital = sqrt(mu/a);
r_sat = [a; 0; 0];      % km (ECI)
v_sat = [0; v_orbital*cos(incl); v_orbital*sin(incl)];  % km/s (ECI)

%% Zaman Parametreleri
% Test tarihi: 1 Ocak 2026, 12:00 UT
year = 2026; month = 1; day = 1;
hour = 12; minute = 0; second = 0;

JD = julian_date(year, month, day, hour, minute, second);
t = 0;  % Epoch'tan beri saniye

fprintf('========================================\n');
fprintf('   ÇEVRE MODELİ TESTİ\n');
fprintf('========================================\n');
fprintf('Tarih: %d/%d/%d %02d:%02d:%02d UT\n', day, month, year, hour, minute, second);
fprintf('Julian Date: %.6f\n', JD);
fprintf('Uydu yüksekliği: %.1f km\n', h_altitude);
fprintf('========================================\n\n');

%% TEST 1: Manyetik Alan (Dipol Model)
fprintf('TEST 1: MANYETİK ALAN (DİPOL MODEL)\n');
fprintf('-----------------------------------\n');

[B_dipole, B_mag_dipole] = magnetic_field_dipole(r_sat, t);

fprintf('Pozisyon (ECI): [%.2f, %.2f, %.2f] km\n', r_sat(1), r_sat(2), r_sat(3));
fprintf('B vektörü (ECI): [%.2f, %.2f, %.2f] nT\n', B_dipole(1), B_dipole(2), B_dipole(3));
fprintf('|B| = %.2f nT (beklenen: ~30,000-50,000 nT @ 1000 km)\n\n', B_mag_dipole);

%% TEST 2: Manyetik Alan (IGRF Model)
fprintf('TEST 2: MANYETİK ALAN (IGRF MODEL)\n');
fprintf('----------------------------------\n');

try
    [B_IGRF, B_NED] = magnetic_field_IGRF(r_sat, t, 6);
    fprintf('B vektörü (ECI): [%.2f, %.2f, %.2f] nT\n', B_IGRF(1), B_IGRF(2), B_IGRF(3));
    fprintf('B vektörü (NED): [%.2f, %.2f, %.2f] nT\n', B_NED(1), B_NED(2), B_NED(3));
    fprintf('|B| (IGRF) = %.2f nT\n', norm(B_IGRF));
    fprintf('Dipol vs IGRF farkı: %.2f nT (%.2f%%)\n\n', ...
            abs(B_mag_dipole - norm(B_IGRF)), ...
            100*abs(B_mag_dipole - norm(B_IGRF))/B_mag_dipole);
catch ME
    fprintf('IGRF hesaplama hatası: %s\n', ME.message);
    fprintf('Dipol model ile devam ediliyor.\n\n');
end

%% TEST 3: Güneş Pozisyonu
fprintf('TEST 3: GÜNEŞ POZİSYONU\n');
fprintf('-----------------------\n');

[r_sun, lambda, epsilon] = sun_position(JD);

fprintf('Güneş pozisyonu (ECI): [%.2e, %.2e, %.2e] km\n', r_sun(1), r_sun(2), r_sun(3));
fprintf('Güneş mesafesi: %.6f AU\n', norm(r_sun) / 149597870.7);
fprintf('Ekliptik boylam: %.2f derece\n', lambda * 180/pi);
fprintf('Ekliptik eğiklik: %.4f derece\n\n', epsilon * 180/pi);

%% TEST 4: Gölge (Eclipse) Hesaplama
fprintf('TEST 4: GÖLGE (ECLIPSE) HESAPLAMA\n');
fprintf('---------------------------------\n');

% Farklı yörünge konumlarında test
test_angles = [0, 90, 180, 270];  % True anomaly (derece)

for i = 1:length(test_angles)
    nu_deg = test_angles(i);
    nu_rad = nu_deg * pi/180;
    
    % Test pozisyonu (dairesel yörünge, equatorial)
    r_test = a * [cos(nu_rad); sin(nu_rad); 0];
    
    [shadow_factor, in_shadow] = shadow_function(r_test, r_sun);
    
    if in_shadow
        status = 'UMBRA (Tam Gölge)';
    elseif shadow_factor < 1.0
        status = sprintf('PENUMBRA (%.1f%%)', shadow_factor*100);
    else
        status = 'GÜNEŞTE';
    end
    
    fprintf('  ν = %3d°: %s\n', nu_deg, status);
end
fprintf('\n');

%% TEST 5: Güneş Radyasyon Basıncı
fprintf('TEST 5: GÜNEŞ RADYASYON BASINCI (SRP)\n');
fprintf('-------------------------------------\n');

[F_srp, a_srp, T_srp] = solar_radiation_pressure(r_sat, r_sun, A_cross, Cr, m);

fprintf('SRP kuvveti: [%.4e, %.4e, %.4e] N\n', F_srp(1), F_srp(2), F_srp(3));
fprintf('SRP ivmesi: [%.4e, %.4e, %.4e] m/s²\n', a_srp(1), a_srp(2), a_srp(3));
fprintf('|F_srp| = %.4e N\n', norm(F_srp));
fprintf('|a_srp| = %.4e m/s² (beklenen: ~10⁻⁸ m/s²)\n\n', norm(a_srp));

%% TEST 6: Atmosferik Sürüklenme
fprintf('TEST 6: ATMOSFERİK SÜRÜKLENME (DRAG)\n');
fprintf('------------------------------------\n');

[a_drag, rho] = atmospheric_drag(r_sat, v_sat, Cd, A_cross, m);

fprintf('Atmosfer yoğunluğu @ %.0f km: %.4e kg/m³\n', h_altitude, rho);
fprintf('Drag ivmesi: [%.4e, %.4e, %.4e] km/s²\n', a_drag(1), a_drag(2), a_drag(3));
fprintf('|a_drag| = %.4e km/s² = %.4e m/s²\n\n', norm(a_drag), norm(a_drag)*1000);

% Farklı yüksekliklerde drag karşılaştırması
fprintf('Yüksekliğe Göre Atmosfer Yoğunluğu:\n');
test_altitudes = [200, 400, 600, 800, 1000];
for i = 1:length(test_altitudes)
    h_test = test_altitudes(i);
    r_test = [(R_earth + h_test); 0; 0];
    [~, rho_test] = atmospheric_drag(r_test, v_sat, Cd, A_cross, m);
    fprintf('  h = %4d km: rho = %.4e kg/m³\n', h_test, rho_test);
end
fprintf('\n');

%% TEST 7: Kütleçekimsel Gradyent Torku
fprintf('TEST 7: KÜTLEÇEKİMSEL GRADYENT TORKU\n');
fprintf('------------------------------------\n');

% Birim dönüşüm matrisi (body = ECI varsayımı)
A_bn = eye(3);

T_gg = gravity_gradient_torque(r_sat, I, A_bn);

fprintf('Atalet momenti: Ix=%.1f, Iy=%.1f, Iz=%.1f kg·m²\n', Ix, Iy, Iz);
fprintf('Gravity gradient torku: [%.4e, %.4e, %.4e] N·m\n', T_gg(1), T_gg(2), T_gg(3));
fprintf('|T_gg| = %.4e N·m (beklenen: ~10⁻⁶ N·m)\n\n', norm(T_gg));

%% Özet Tablosu
fprintf('========================================\n');
fprintf('   ÇEVRE ETKİLERİ ÖZET TABLOSU\n');
fprintf('========================================\n');
fprintf('Etki                    | Büyüklük\n');
fprintf('------------------------|-----------------\n');
fprintf('Manyetik alan           | %.0f nT\n', B_mag_dipole);
fprintf('SRP ivmesi              | %.2e m/s²\n', norm(a_srp));
fprintf('Drag ivmesi             | %.2e m/s²\n', norm(a_drag)*1000);
fprintf('Gravity gradient torku  | %.2e N·m\n', norm(T_gg));
fprintf('========================================\n');

%% Grafikler
fprintf('\nGRAFİKLER OLUŞTURULUYOR...\n');

figure('Name', 'Çevre Modeli Sonuçları', 'Position', [100 100 1200 800]);

% 1. Manyetik Alan vs Yükseklik
subplot(2,2,1);
altitudes = linspace(200, 2000, 50);
B_values = zeros(size(altitudes));
for i = 1:length(altitudes)
    r_test = [(R_earth + altitudes(i)); 0; 0];
    [~, B_values(i)] = magnetic_field_dipole(r_test, t);
end
plot(altitudes, B_values/1000, 'b-', 'LineWidth', 2);
xlabel('Yükseklik (km)');
ylabel('|B| (μT)');
title('Manyetik Alan Şiddeti vs Yükseklik');
grid on;

% 2. Atmosfer Yoğunluğu vs Yükseklik
subplot(2,2,2);
altitudes_atm = linspace(100, 1000, 50);
rho_values = zeros(size(altitudes_atm));
for i = 1:length(altitudes_atm)
    r_test = [(R_earth + altitudes_atm(i)); 0; 0];
    [~, rho_values(i)] = atmospheric_drag(r_test, v_sat, Cd, A_cross, m);
end
semilogy(altitudes_atm, rho_values, 'r-', 'LineWidth', 2);
xlabel('Yükseklik (km)');
ylabel('Yoğunluk (kg/m³)');
title('Atmosfer Yoğunluğu vs Yükseklik');
grid on;

% 3. Güneş-Dünya-Uydu Geometrisi (XY düzlemi)
subplot(2,2,3);
theta_orbit = linspace(0, 2*pi, 100);
x_orbit = a * cos(theta_orbit);
y_orbit = a * sin(theta_orbit);
plot(x_orbit, y_orbit, 'b-', 'LineWidth', 1.5);
hold on;
% Dünya
theta_earth = linspace(0, 2*pi, 100);
plot(R_earth*cos(theta_earth), R_earth*sin(theta_earth), 'g-', 'LineWidth', 2);
% Uydu
plot(r_sat(1), r_sat(2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% Güneş yönü
sun_dir = r_sun / norm(r_sun);
quiver(0, 0, sun_dir(1)*R_earth*2, sun_dir(2)*R_earth*2, 'y', 'LineWidth', 2);
xlabel('X (km)'); ylabel('Y (km)');
title('Yörünge ve Güneş Yönü (XY Düzlemi)');
axis equal; grid on;
legend('Yörünge', 'Dünya', 'Uydu', 'Güneş Yönü', 'Location', 'best');

% 4. Çevre Etkileri Karşılaştırması
subplot(2,2,4);
effects = categorical({'SRP', 'Drag', 'J2'});
values = [norm(a_srp), norm(a_drag)*1000, 1e-5];  % m/s² (J2 yaklaşık)
bar(effects, log10(values));
ylabel('log_{10}(İvme) [m/s²]');
title('Pertürbasyon İvmeleri Karşılaştırması');
grid on;

fprintf('Tüm testler tamamlandı.\n');
fprintf('========================================\n');
