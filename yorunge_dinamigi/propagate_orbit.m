% Orbital Mechanics for Engineering Students kitabından uyarlanmıştır.
% ODE entegrasyonu ile yörünge propagasyonu

% Bu fonksiyon, ODE çözücü kullanarak yörüngeyi propagate eder.
% J2 pertürbasyonu dahil veya hariç çalışabilir.

% r0        - başlangıç pozisyon vektörü (km) [1x3]
% v0        - başlangıç hız vektörü (km/s) [1x3]
% tspan     - zaman aralığı [t_start, t_end] veya zaman vektörü (s)
% options   - ODE çözücü seçenekleri (opsiyonel)
% use_J2    - J2 pertürbasyonu kullan (true/false, varsayılan: false)
% t_out     - çıkış zaman vektörü (s)
% r_out     - çıkış pozisyon matrisi (km) [Nx3]
% v_out     - çıkış hız matrisi (km/s) [Nx3]

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: 
%   two_body_EOM, two_body_J2_EOM, J2_perturbation

function [t_out, r_out, v_out] = propagate_orbit(r0, v0, tspan, use_J2)

% Varsayılan: J2 kullanma
if nargin < 4
    use_J2 = false;
end

% Başlangıç state vektörü (sütun vektörü olmalı)
if isrow(r0)
    r0 = r0';
end
if isrow(v0)
    v0 = v0';
end

state0 = [r0; v0];

% ODE çözücü seçenekleri (yüksek hassasiyet)
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

% Hareket denklemini seç
if use_J2
    EOM = @two_body_J2_EOM;
else
    EOM = @two_body_EOM;
end

% ODE çözümü
[t_out, state_out] = ode113(EOM, tspan, state0, options);

% Pozisyon ve hız matrislerini ayır
r_out = state_out(:, 1:3);
v_out = state_out(:, 4:6);

end
