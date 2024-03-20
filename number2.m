function Tp = nilaitp()
    Cp = 0.005;
    D = 100 * 10^6;
    teta = 0.2;
    b = 69.4;
    alpa = 3;
    k = 0.12 * 365;
    Cb = 0.0315;
    betaa = 81;
    A = 3200;
    Ch = 0.002;
    CO = 5000;

    f = @(Tp) Cp / (1 + b) * (alpa * (1 + b * exp(-k * Tp)) - b * k * exp(-k * Tp))...
        + Cb / (A * betaa) * (alpa * (1 + b * exp(-k * Tp)) * (exp(betaa * Tp) - 1)...
        - b * k * exp(-k * Tp) * (exp(betaa * Tp) - 1) + betaa * (1 + b * exp(-k * Tp)) * exp(betaa * Tp));

    Tp0 = 0;
    options = optimset('Display', 'iter', 'TolFun', 1e-8, 'TolX', 1e-8);
    Tp = fsolve(f, Tp0, options);
end
