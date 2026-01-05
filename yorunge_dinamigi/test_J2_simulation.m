% =========================================================================
% LEO Uydu Yörünge Simülasyonu - J2 Pertürbasyonu Dahil
% 60 kg uydu, 1000 km yükseklik
% =========================================================================
% Bu script, J2 pertürbasyonunun yörünge üzerindeki etkilerini gösterir.
% İki cisim modeli ile J2 dahil modeli karşılaştırır.
%
% Kullanılan fonksiyonlar:
%   - sv_from_coe.m        : Orbital elementlerden state vektör hesapla
%   - coe_from_sv.m        : State vektörden orbital elementler hesapla
%   - propagate_orbit.m    : ODE tabanlı yörünge propagatörü
%   - two_body_EOM.m       : İki cisim hareket denklemi
%   - two_body_J2_EOM.m    : J2 dahil hareket denklemi
%   - J2_perturbation.m    : J2 ivme hesaplama
%   - J2_secular_rates.m   : Analitik RAAN/omega drift hızları
% =========================================================================

clear; clc; close all;

%% Sabitler
global mu
mu = 398600;            % Dünya kütleçekimsel parametresi (km^3/s^2)
R_earth = 6378.137;     % Dünya yarıçapı (km)
J2 = 1.08263e-3;        % J2 katsayısı

%% Yörünge Parametreleri (1000 km LEO, Güneş-senkron)
h_altitude = 1000;                      % Yörünge yüksekliği (km)
a = R_earth + h_altitude;               % Yarı-büyük eksen (km)
e = 0.0001;                             % Küçük eksantriklik (tam 0 sorun çıkarabilir)
incl = 98 * pi/180;                     % Eğilim (rad) - Güneş-senkron için tipik
RA = 0 * pi/180;                        % RAAN (rad)
w = 0 * pi/180;                         % Argument of perigee (rad)
TA = 0 * pi/180;                        % True anomaly - başlangıç (rad)

% Dairesel yörünge için açısal momentum hesapla
h = sqrt(mu * a * (1 - e^2));           % km^2/s

% Yörünge periyodu
T = 2*pi*sqrt(a^3/mu);                  % saniye
T_min = T / 60;                         % dakika
T_hour = T / 3600;                      % saat

% Orbital hız (dairesel)
v_orbital = sqrt(mu/a);                 % km/s

%% Yörünge Bilgilerini Ekrana Yazdır
fprintf('========================================\n');
fprintf('   LEO UYDU YÖRÜNGE PARAMETRELERİ\n');
fprintf('========================================\n');
fprintf('Yörünge yüksekliği  : %.1f km\n', h_altitude);
fprintf('Yarı-büyük eksen    : %.3f km\n', a);
fprintf('Eksantriklik        : %.6f\n', e);
fprintf('Eğilim              : %.2f°\n', incl*180/pi);
fprintf('RAAN                : %.2f°\n', RA*180/pi);
fprintf('Yörünge periyodu    : %.2f dakika (%.4f saat)\n', T_min, T_hour);
fprintf('Orbital hız         : %.4f km/s\n', v_orbital);
fprintf('========================================\n\n');

%% Analitik J2 Etkileri (Seküler Drift Hızları)
fprintf('ANALİTİK J2 SEKülER ETKİLERİ\n');
fprintf('----------------------------\n');

[RAAN_dot, w_dot] = J2_secular_rates(a, e, incl);

% Birim dönüşümleri
RAAN_dot_deg_day = RAAN_dot * 180/pi * 86400;   % derece/gün
w_dot_deg_day = w_dot * 180/pi * 86400;         % derece/gün

fprintf('RAAN drift hızı     : %.6f °/gün\n', RAAN_dot_deg_day);
fprintf('Omega drift hızı    : %.6f °/gün\n', w_dot_deg_day);

% Güneş-senkron kontrolü (yaklaşık 0.9856°/gün olmalı)
fprintf('\nGüneş-senkron için gereken RAAN drift: ~0.9856 °/gün\n');
if abs(RAAN_dot_deg_day - 0.9856) < 0.1
    fprintf('✓ Bu yörünge güneş-senkron yakınında!\n');
else
    fprintf('✗ Güneş-senkron için eğilimi ayarla.\n');
end
fprintf('\n');

%% Başlangıç State Vektörü
fprintf('BAŞLANGIÇ STATE VEKTÖRÜ HESAPLAMA\n');
fprintf('---------------------------------\n');

% Orbital elementler vektörü [h e RA incl w TA]
coe = [h, e, RA, incl, w, TA];

% State vektörü hesapla
[r0, v0] = sv_from_coe(coe);

fprintf('r0 = [%.4f, %.4f, %.4f] km\n', r0(1), r0(2), r0(3));
fprintf('v0 = [%.4f, %.4f, %.4f] km/s\n', v0(1), v0(2), v0(3));
fprintf('|r0| = %.4f km\n', norm(r0));
fprintf('|v0| = %.4f km/s\n\n', norm(v0));

%% Simülasyon Ayarları
fprintf('SİMÜLASYON AYARLARI\n');
fprintf('-------------------\n');

num_orbits = 15;                        % Simüle edilecek yörünge sayısı
t_total = num_orbits * T;               % Toplam süre (s)
num_points = num_orbits * 100;          % Toplam nokta sayısı
tspan = linspace(0, t_total, num_points);

fprintf('Simülasyon süresi   : %.2f saat (%.1f yörünge)\n', t_total/3600, num_orbits);
fprintf('Zaman adımı sayısı  : %d\n\n', num_points);

%% Yörünge Propagasyonu - İki Cisim (J2 olmadan)
fprintf('YÖRÜNGE PROPAGASYONU\n');
fprintf('--------------------\n');
fprintf('1. İki cisim modeli (J2 yok)...\n');
tic;
[t_2body, r_2body, v_2body] = propagate_orbit(r0, v0, tspan, false);
time_2body = toc;
fprintf('   Tamamlandı (%.2f s)\n', time_2body);

%% Yörünge Propagasyonu - J2 Dahil
fprintf('2. J2 pertürbasyonlu model...\n');
tic;
[t_J2, r_J2, v_J2] = propagate_orbit(r0, v0, tspan, true);
time_J2 = toc;
fprintf('   Tamamlandı (%.2f s)\n\n', time_J2);

%% Orbital Elementlerin Zaman Evrimi
fprintf('ORBİTAL ELEMENT EVRİMİ HESAPLANIYOR...\n');

% Her 10 noktada bir orbital elementleri hesapla (hız için)
sample_rate = 10;
idx_samples = 1:sample_rate:length(t_J2);
num_samples = length(idx_samples);

% Sonuç dizileri
coe_2body = zeros(num_samples, 7);
coe_J2 = zeros(num_samples, 7);

for i = 1:num_samples
    idx = idx_samples(i);
    coe_2body(i,:) = coe_from_sv(r_2body(idx,:), v_2body(idx,:));
    coe_J2(i,:) = coe_from_sv(r_J2(idx,:), v_J2(idx,:));
end

t_samples = t_J2(idx_samples);

fprintf('Tamamlandı.\n\n');

%% Pozisyon Farkı Hesaplama
dr = r_J2 - r_2body;
dr_mag = sqrt(dr(:,1).^2 + dr(:,2).^2 + dr(:,3).^2);

fprintf('POZİSYON FARKI (J2 vs İki-Cisim)\n');
fprintf('--------------------------------\n');
fprintf('Maksimum fark: %.4f km\n', max(dr_mag));
fprintf('Son andaki fark: %.4f km\n', dr_mag(end));
fprintf('%.1f yörünge sonra fark: %.2f km\n\n', num_orbits, dr_mag(end));

%% Grafikler
fprintf('GRAFİKLER OLUŞTURULUYOR...\n');

% ===================== Şekil 1: 3D Yörünge Karşılaştırma =====================
figure('Name', 'J2 Etkisi - 3D Yörünge', 'Position', [50 50 1000 800]);

% Dünya'yı çiz
[X_earth, Y_earth, Z_earth] = sphere(50);
surf(X_earth*R_earth, Y_earth*R_earth, Z_earth*R_earth, ...
     'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
hold on;

% İki cisim yörüngesi
plot3(r_2body(:,1), r_2body(:,2), r_2body(:,3), 'b-', 'LineWidth', 1);

% J2 yörüngesi
plot3(r_J2(:,1), r_J2(:,2), r_J2(:,3), 'r-', 'LineWidth', 1);

% Başlangıç noktası
plot3(r0(1), r0(2), r0(3), 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');

% ECI eksenleri
axis_len = R_earth * 1.3;
quiver3(0,0,0, axis_len, 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3);
quiver3(0,0,0, 0, axis_len, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3);
quiver3(0,0,0, 0, 0, axis_len, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3);
text(axis_len*1.1, 0, 0, 'X', 'FontSize', 12);
text(0, axis_len*1.1, 0, 'Y', 'FontSize', 12);
text(0, 0, axis_len*1.1, 'Z', 'FontSize', 12);

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title(sprintf('LEO Yörüngesi: İki-Cisim vs J2 (%d yörünge)', num_orbits));
legend('Dünya', 'İki-Cisim', 'J2 Dahil', 'Başlangıç', 'Location', 'best');
axis equal; grid on;
view(45, 30);

% ===================== Şekil 2: Orbital Element Değişimleri =====================
figure('Name', 'J2 Etkisi - Orbital Elementler', 'Position', [100 100 1200 800]);

% RAAN değişimi
subplot(2,2,1);
plot(t_samples/3600, coe_2body(:,3)*180/pi, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_samples/3600, coe_J2(:,3)*180/pi, 'r-', 'LineWidth', 1.5);
% Analitik tahmin
RAAN_analytic = RA*180/pi + RAAN_dot_deg_day * t_samples/86400;
plot(t_samples/3600, RAAN_analytic, 'g--', 'LineWidth', 1.5);
xlabel('Zaman (saat)');
ylabel('RAAN (°)');
title('RAAN (Ω) Değişimi');
legend('İki-Cisim', 'J2 Dahil', 'Analitik J2', 'Location', 'best');
grid on;

% Argument of Perigee değişimi
subplot(2,2,2);
plot(t_samples/3600, coe_2body(:,5)*180/pi, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_samples/3600, coe_J2(:,5)*180/pi, 'r-', 'LineWidth', 1.5);
% Analitik tahmin
w_analytic = w*180/pi + w_dot_deg_day * t_samples/86400;
plot(t_samples/3600, w_analytic, 'g--', 'LineWidth', 1.5);
xlabel('Zaman (saat)');
ylabel('ω (°)');
title('Argument of Perigee (ω) Değişimi');
legend('İki-Cisim', 'J2 Dahil', 'Analitik J2', 'Location', 'best');
grid on;

% Eğilim değişimi
subplot(2,2,3);
plot(t_samples/3600, coe_2body(:,4)*180/pi, 'b-', 'LineWidth', 1.5);
hold on;
plot(t_samples/3600, coe_J2(:,4)*180/pi, 'r-', 'LineWidth', 1.5);
xlabel('Zaman (saat)');
ylabel('Eğilim (°)');
title('Eğilim (i) Değişimi');
legend('İki-Cisim', 'J2 Dahil', 'Location', 'best');
grid on;

% Pozisyon farkı
subplot(2,2,4);
plot(t_J2/3600, dr_mag, 'k-', 'LineWidth', 1.5);
xlabel('Zaman (saat)');
ylabel('Pozisyon Farkı (km)');
title('|r_{J2} - r_{2body}| Zamanla Değişim');
grid on;

% ===================== Şekil 3: Ground Track =====================
figure('Name', 'J2 Etkisi - Ground Track', 'Position', [150 150 1200 500]);

% Basit enlem/boylam hesaplama (ECEF dönüşümü olmadan yaklaşık)
% Dünya'nın dönüş hızı
omega_earth = 7.2921159e-5;  % rad/s

% J2 yörüngesi için ground track
lon_J2 = zeros(length(t_J2), 1);
lat_J2 = zeros(length(t_J2), 1);

for i = 1:length(t_J2)
    % ECI'dan yaklaşık ECEF'e (Dünya dönüşü)
    theta = omega_earth * t_J2(i);
    x_ecef = r_J2(i,1)*cos(theta) + r_J2(i,2)*sin(theta);
    y_ecef = -r_J2(i,1)*sin(theta) + r_J2(i,2)*cos(theta);
    z_ecef = r_J2(i,3);
    
    % Enlem ve boylam
    r_xy = sqrt(x_ecef^2 + y_ecef^2);
    lat_J2(i) = atan2(z_ecef, r_xy) * 180/pi;
    lon_J2(i) = atan2(y_ecef, x_ecef) * 180/pi;
end

% Ground track çiz
plot(lon_J2, lat_J2, 'r.', 'MarkerSize', 1);
hold on;

% Dünya haritası (basit sınırlar)
plot([-180 180], [0 0], 'k-', 'LineWidth', 0.5);  % Ekvator
plot([0 0], [-90 90], 'k-', 'LineWidth', 0.5);    % Başlangıç meridyeni

xlabel('Boylam (°)');
ylabel('Enlem (°)');
title(sprintf('Ground Track - %d Yörünge (J2 Dahil)', num_orbits));
xlim([-180 180]);
ylim([-90 90]);
grid on;
axis equal;

%% Özet
fprintf('\n========================================\n');
fprintf('   SİMÜLASYON TAMAMLANDI\n');
fprintf('========================================\n');
fprintf('J2 pertürbasyonunun etkileri:\n');
fprintf('  • RAAN kayması: %.4f °/gün\n', RAAN_dot_deg_day);
fprintf('  • ω kayması: %.4f °/gün\n', w_dot_deg_day);
fprintf('  • %d yörünge sonra pozisyon farkı: %.2f km\n', num_orbits, dr_mag(end));
fprintf('========================================\n');
fprintf('Sonraki adım: Simulink modeli oluşturmak.\n');
fprintf('========================================\n');
