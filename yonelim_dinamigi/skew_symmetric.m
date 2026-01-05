% Spacecraft Attitude Determination and Control (Wertz) kitabından alınmıştır.
% Denklem 16-6b: Skew-Symmetric Matrix

% Bu fonksiyon, bir 3x1 vektörden skew-symmetric (çapraz çarpım) matrisi oluşturur.

% v         - 3x1 vektör [v1, v2, v3]'
% S         - 3x3 skew-symmetric matris

% Skew-Symmetric Matris:
%       [  0   -v3   v2  ]
%   S = [  v3   0   -v1  ]
%       [ -v2   v1   0   ]

% Özellik: S * u = v × u (çapraz çarpım)

% kullanıcı tarafından tanımlanmış ihtiyaç olan fonksiyonlar: yok

function S = skew_symmetric(v)

% Vektörü sütun vektörüne çevir
if isrow(v)
    v = v';
end

% Skew-symmetric matris
S = [  0     -v(3)   v(2);
       v(3)   0     -v(1);
      -v(2)   v(1)   0   ];

end
