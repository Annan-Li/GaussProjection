FE = 500000;
FN = 0;
k0 = 1;
L0 = 111 / 180 * pi;%中央子午线经度（弧度）
% a = 6378137;%长半轴 国家2000坐标系
% f = 1 / 298.257222101;%扁率

a = 6378245;%长半轴 北京54坐标系
f = 1 / 298.3;%扁率

X = 2505530.1986;
Y = 710280.0416;

e = (2 * f - f * f) ^ 0.5;
ep = ((1 / (1 - f) ^ 2) - 1) ^ 0.5;%参考椭球第二偏心率

ef = f / (2 - f);

Mf = (X - FN) / k0;

B1 = Mf / (a * (1 - e * e / 4 - 3 * e * e * e * e / 64 - 5 * e * e * e * e * e * e / 256));

Bf = B1 + (3 * ef / 2 - 27 * ef * ef * ef / 32) * sin(2 * B1) + (21 * ef * ef / 16 - 55 * ef * ef * ef * ef / 32) * sin(4 * B1) + 151 * ef * ef * ef / 96 * sin(6 * B1) + 1097 * ef * ef * ef * ef / 512 * sin(8 * B1);

Nf = a / sqrt(1 - e * e * sin(Bf) * sin(Bf));
Rf = a * (1 - e * e) / (1 - e * e * sin(Bf) * sin(Bf)) ^ 1.5;
Tf = tan(Bf) * tan(Bf);
Cf = ep * ep * cos(Bf) * cos(Bf);
D = (Y - FE) / (Nf * k0);

L = L0 +(1 / cos(Bf)) * (D - (1 + 2 * Tf + Cf) * D * D * D / 6 +(5 - 2 * Cf + 28 * Tf - 3 * Cf * Cf + 8 * ep * ep + 24 * Tf * Tf) * D * D * D * D * D / 120)
B = Bf - Nf * tan(Bf) / Rf * (D * D / 2 - (5 + 3 * Tf + 10 * Cf - 4 * Cf * Cf - 9 * ep * ep) * D * D * D * D / 24 + (61 + 90 * Tf + 298 * Cf + 45 * Tf * Tf - 252 * ep * ep - 3 * Cf * Cf) * D * D * D * D * D * D / 720)



