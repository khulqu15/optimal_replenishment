function optimization1()
    A = 3200;
    b = 69.4;
    k = 0.12 * 365; % Mengubah basis waktu menjadi tahun
    cp = 0.005;
    cb = 0.0315;    
    ch = 0.002;
    co = 5000;
    D = 100 * 10^6; % Laju permintaan konstan (g/tahun)
    theta = 0.2;
    beta = 81; 
    alpha = 3;

    TD = 0.2146;
    TP = 0.0719;

    TUC = evaluate_expression(cp, D, theta, b, TD, alpha, TP, k, cb, beta, A, ch, co);
    
    fprintf('TP (Periode Breeding): %.4f tahun\n', TP);
    fprintf('TD (Periode Konsumsi): %.4f tahun\n', TD);
    fprintf('Total Unit Cost (TUC): %.2f€\n', TUC);

    function z = evaluate_expression(Cp, D, teta, b, Td, alpa, Tp, k, Cb, betaa, A, Ch, CO)
        Q = D / (teta * (1 + b)) * (exp(teta * Td) - 1) * exp(alpa * Tp) * (1 + b * exp(-k * Tp));
        PC = (Cp * D) / (teta * (   1 + b)) * (exp(teta * Td) - 1) * exp(alpa * Tp) * (1 + b * exp(-k * Tp));
        BC = (Cb * D) / (teta * betaa * A) * (exp(teta * Td) - 1) * exp(alpa * Tp) * (1 + b * exp(-k * Tp)) * (exp(betaa * Tp) - 1);
        HC = Ch * (D / teta) * ((1 / teta) * exp(teta * Td) - (1 / teta) - Td);

        OC = CO;
        z = (PC + BC + HC + OC) / Td;
    end
    
    function F = systemOfEquations(TP, cp, D, theta, b, alpha, k, cb, beta, A, ch, co)
    F = cp / (1 + b) * (alpha * (1 + b * exp(-k * TP)) - b * k * exp(-k * TP)) + ...
        cb / (A * beta) * (alpha * (1 + b * exp(-k * TP)) * (exp(beta * TP) -xcorr 1)) - ...
        (b * k * exp(-k * TP) * (exp(beta * TP) - 1) + beta * (1 + b * exp(-k * TP)) * exp(beta * TP));
    end
end
