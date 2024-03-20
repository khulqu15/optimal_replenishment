Cp = 0.005;
D = 100000000;
teta = 0.2;
b = 69.4;
alpa = 3;
k = 43.8;
Cb = 0.0315;
betaa = 81;
A = 3200;
Ch = 0.002;
CO = 5000;

f = @(Tp) (Cp / (1 + b) * (alpa * (1 + b * exp(-k * Tp)) - b * k * exp(-k * Tp))) ...
    + (Cb / (betaa * A) * (alpa * (1 + b * exp(-k * Tp)) * (exp(betaa * Tp) - 1))) ...
    - (b * k * exp(-k * Tp) * (exp(betaa * Tp) - 1) + betaa * (1 + b * exp(-k * Tp)) * exp(betaa * Tp));

Tp0 = linspace(0.5, 0.8, 100);
fValues = arrayfun(f, TpRange);

options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8);

optimizedTp = zeros(size(TpRange));
optimizedFValues = zeros(size(TpRange));

for i = 1:length(TpRange)
    [optimizedTp(i), fval] = fsolve(f, TpRange(i), options);
    optimizedFValues(i) = fval;
    
    fprintf('Tp = %.4f, f(Tp) = %.4f, Tp0 = %.4f\n', optimizedTp(i), optimizedFValues(i), TpRange(i));
end

figure;
plot(optimizedTp, TpRange, 'ro');
xlabel('Tp');
ylabel('Tp0');
title('Hasil Optimasi');
grid on;
