% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Chapter 16 - Attitude Dynamics
% Denklem 16-1: Quaternion Propagation (Discrete)

% Bu fonksiyon, quaternion'u zaman adımında ilerletir.
% Sabit açısal hız varsayımı ile tam çözüm kullanır.

% q         - mevcut quaternion [4x1]
% omega     - açısal hız vektörü body frame'de (rad/s) [3x1]
% dt        - zaman adımı (s)
% q_new     - yeni quaternion [4x1]

% Wertz Eq. 16-1:
%   q(t+Δt) = [cos(Δθ/2)*I + sin(Δθ/2)/θ * Ω] * q(t)
%
% Burada Δθ = |ω| * Δt

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function q_new = quat_propagate(q, omega, dt)

% Vektörleri sütun vektörüne çevir
if isrow(q)
    q = q';
end
if isrow(omega)
    omega = omega';
end

% Açısal hız büyüklüğü
omega_mag = norm(omega);

% Eğer açısal hız çok küçükse, quaternion değişmez
if omega_mag < 1e-12
    q_new = q;
    return;
end

% Dönme açısı
d_theta = omega_mag * dt;

% Dönme ekseni birim vektörü
e = omega / omega_mag;

% Delta quaternion (küçük dönüşü temsil eder)
% q_delta = [e * sin(Δθ/2); cos(Δθ/2)]
sin_half = sin(d_theta / 2);
cos_half = cos(d_theta / 2);

q_delta = [e(1) * sin_half;
           e(2) * sin_half;
           e(3) * sin_half;
           cos_half];

% Quaternion çarpımı: q_new = q_delta ⊗ q
% Bu, önce q sonra q_delta dönüşümünü verir
q_new = quat_mult_local(q_delta, q);

% Normalize et
q_new = q_new / norm(q_new);

end


function q = quat_mult_local(q1, q2)
% Yerel quaternion çarpımı (bağımlılık önlemek için)
a1 = q1(1); b1 = q1(2); c1 = q1(3); d1 = q1(4);
a2 = q2(1); b2 = q2(2); c2 = q2(3); d2 = q2(4);

q = [d1*a2 + a1*d2 + b1*c2 - c1*b2;
     d1*b2 - a1*c2 + b1*d2 + c1*a2;
     d1*c2 + a1*b2 - b1*a2 + c1*d2;
     d1*d2 - a1*a2 - b1*b2 - c1*c2];
end
