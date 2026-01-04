% =========================================================================
% LEO Uydu Yörünge Simülasyonu - Test Script
% 60 kg uydu, 1000 km yükseklik
% =========================================================================
% Bu script, tüm yörünge fonksiyonlarını test eder ve basit bir simülasyon
% çalıştırır.
%
% Kullanılan fonksiyonlar:
%   - sv_from_coe.m    : Orbital elementlerden state vektör hesapla
%   - rv_from_r0v0.m   : State vektör propagasyonu
%   - coe_from_sv.m    : State vektörden orbital elementler hesapla
%   - kepler_U.m       : Universal Kepler denklemi çözücü
%   - f_and_g.m        : Lagrange f ve g katsayıları
%   - fDot_and_gDot.m  : Lagrange katsayılarının türevleri
%   - stumpC.m         : Stumpff C fonksiyonu
%   - stumpS.m         : Stumpff S fonksiyonu
% =========================================================================

clear; clc; close all;

%% Sabitler
global mu
mu = 398600;            % Dünya kütleçekimsel parametresi (km^3/s^2)
R_earth = 6378.137;     % Dünya yarıçapı (km)

%% Yörünge Parametreleri (1000 km LEO, Güneş-senkron)
h_altitude = 1000;                      % Yörünge yüksekliği (km)
a = R_earth + h_altitude;               % Yarı-büyük eksen (km)
e = 0;                                  % Eksantriklik (dairesel yörünge)
incl = 98 * pi/180;                     % Eğilim (rad) - Güneş-senkron için tipik
RA = 0 * pi/180;                        % RAAN (rad)
w = 0 * pi/180;                         % Argument of perigee (rad)
TA = 0 * pi/180;                        % True anomaly - başlangıç (rad)

% Dairesel yörünge için açısal momentum hesapla
h = sqrt(mu * a * (1 - e^2));           % km^2/s

% Yörünge periyodu
T = 2*pi*sqrt(a^3/mu);                  % saniye
T_min = T / 60;                         % dakika

% Orbital hız (dairesel)
v_orbital = sqrt(mu/a);                 % km/s

%% Yörünge Bilgilerini Ekrana Yazdır
fprintf('========================================\n');
fprintf('   LEO UYDU YÖRÜNGE PARAMETRELERİ\n');
fprintf('========================================\n');
fprintf('Yörünge yüksekliği  : %.1f km\n', h_altitude);
fprintf('Yarı-büyük eksen    : %.3f km\n', a);
fprintf('Eksantriklik        : %.4f\n', e);
fprintf('Eğilim              : %.2f°\n', incl*180/pi);
fprintf('RAAN                : %.2f°\n', RA*180/pi);
fprintf('Yörünge periyodu    : %.2f dakika\n', T_min);
fprintf('Orbital hız         : %.4f km/s\n', v_orbital);
fprintf('Açısal momentum     : %.4f km^2/s\n', h);
fprintf('========================================\n\n');

%% ADIM 1: Orbital Elementlerden State Vektör Hesapla
fprintf('ADIM 1: sv_from_coe() ile State Vektör Hesaplama\n');
fprintf('------------------------------------------------\n');

% Orbital elementler vektörü [h e RA incl w TA]
coe = [h, e, RA, incl, w, TA];

% State vektörü hesapla
[r0, v0] = sv_from_coe(coe);

fprintf('Başlangıç pozisyonu r0 = [%.4f, %.4f, %.4f] km\n', r0(1), r0(2), r0(3));
fprintf('Başlangıç hızı     v0 = [%.4f, %.4f, %.4f] km/s\n', v0(1), v0(2), v0(3));
fprintf('|r0| = %.4f km (beklenen: %.4f km)\n', norm(r0), a);
fprintf('|v0| = %.4f km/s (beklenen: %.4f km/s)\n\n', norm(v0), v_orbital);

%% ADIM 2: Yörünge Propagasyonu (1 tam periyot)
fprintf('ADIM 2: rv_from_r0v0() ile Yörünge Propagasyonu\n');
fprintf('-----------------------------------------------\n');

% Simülasyon ayarları
num_steps = 100;                        % Adım sayısı
t_total = T;                            % Toplam süre (1 periyot)
dt = t_total / num_steps;               % Zaman adımı

% Sonuçları saklamak için diziler
t_vec = zeros(1, num_steps+1);
r_vec = zeros(num_steps+1, 3);
v_vec = zeros(num_steps+1, 3);

% Başlangıç değerleri
t_vec(1) = 0;
r_vec(1,:) = r0;
v_vec(1,:) = v0;

% Propagasyon döngüsü
fprintf('Propagasyon başlıyor (%d adım)...\n', num_steps);
for i = 1:num_steps
    t = i * dt;
    [r_new, v_new] = rv_from_r0v0(r0, v0, t);
    t_vec(i+1) = t;
    r_vec(i+1,:) = r_new;
    v_vec(i+1,:) = v_new;
end
fprintf('Propagasyon tamamlandı.\n\n');

% Son değerleri kontrol et (1 periyot sonra başlangıca dönmeli)
fprintf('1 periyot sonra:\n');
fprintf('Son pozisyon r = [%.4f, %.4f, %.4f] km\n', r_vec(end,1), r_vec(end,2), r_vec(end,3));
fprintf('Başlangıç    r0= [%.4f, %.4f, %.4f] km\n', r0(1), r0(2), r0(3));
fprintf('Fark = %.6f km (sıfıra yakın olmalı)\n\n', norm(r_vec(end,:) - r0));

%% ADIM 3: State Vektörden Orbital Elementler (Doğrulama)
fprintf('ADIM 3: coe_from_sv() ile Orbital Element Doğrulama\n');
fprintf('---------------------------------------------------\n');

% Yarı periyottaki state vektörü al
mid_idx = round(num_steps/2);
r_mid = r_vec(mid_idx,:);
v_mid = v_vec(mid_idx,:);

% Orbital elementleri hesapla
coe_check = coe_from_sv(r_mid, v_mid);

fprintf('Orijinal vs Hesaplanan:\n');
fprintf('  h    : %.4f vs %.4f km^2/s\n', h, coe_check(1));
fprintf('  e    : %.6f vs %.6f\n', e, coe_check(2));
fprintf('  RAAN : %.4f° vs %.4f°\n', RA*180/pi, coe_check(3)*180/pi);
fprintf('  incl : %.4f° vs %.4f°\n', incl*180/pi, coe_check(4)*180/pi);
fprintf('  a    : %.4f vs %.4f km\n\n', a, coe_check(7));

%% ADIM 4: 3D Yörünge Grafiği
fprintf('ADIM 4: Yörünge Görselleştirme\n');
fprintf('-----------------------------\n');

figure('Name', 'LEO Uydu Yörüngesi', 'Position', [100 100 1200 500]);

% 3D Yörünge
subplot(1,2,1);
plot3(r_vec(:,1), r_vec(:,2), r_vec(:,3), 'b-', 'LineWidth', 1.5);
hold on;

% Dünya'yı çiz (basit küre)
[X_earth, Y_earth, Z_earth] = sphere(30);
surf(X_earth*R_earth, Y_earth*R_earth, Z_earth*R_earth, ...
     'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

% Başlangıç noktası
plot3(r0(1), r0(2), r0(3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

% ECI eksenleri
quiver3(0,0,0, R_earth*1.5, 0, 0, 'r', 'LineWidth', 2);
quiver3(0,0,0, 0, R_earth*1.5, 0, 'g', 'LineWidth', 2);
quiver3(0,0,0, 0, 0, R_earth*1.5, 'b', 'LineWidth', 2);

xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
title(sprintf('1000 km LEO Yörüngesi (i=%.0f°)', incl*180/pi));
axis equal; grid on;
legend('Yörünge', 'Dünya', 'Başlangıç', 'Location', 'best');
view(45, 30);

% Yükseklik vs Zaman
subplot(1,2,2);
altitude = sqrt(r_vec(:,1).^2 + r_vec(:,2).^2 + r_vec(:,3).^2) - R_earth;
plot(t_vec/60, altitude, 'b-', 'LineWidth', 1.5);
xlabel('Zaman (dakika)');
ylabel('Yükseklik (km)');
title('Yükseklik vs Zaman');
grid on;
ylim([h_altitude-10, h_altitude+10]);

fprintf('Grafikler oluşturuldu.\n\n');

%% Özet
fprintf('========================================\n');
fprintf('   SİMÜLASYON TAMAMLANDI\n');
fprintf('========================================\n');
fprintf('Tüm fonksiyonlar başarıyla çalıştı!\n');
fprintf('Sonraki adım: J2 pertürbasyonu eklemek.\n');
fprintf('========================================\n');
