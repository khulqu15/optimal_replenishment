function optimized_solution()
    % Definisi parameter masalah
    D = 100 * 10^6;
    cp = 0.005;
    lb = [0.222, 0.0015];
    ub = [1, 2];
    cb = 0.0315;
    ch = 0.002;
    co = 5000;
    theta = 0.2;
    beta = 81; 
    alpha = 3;
    k = 0.12 * 365;
    A = 3200;
    delta = 850;
    TP = 0.0719;
    b = 69.4;
    xi = 0.0015;
    x0 = 0.00;

    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ... % Mengatur algoritma yang digunakan oleh fmincon. 'sqp' berarti Sequential Quadratic Programming, yang efektif untuk non-linear optimization.
        'Display', 'iter'); % Mengatur tingkat tampilan output. 'iter' akan menampilkan informasi setiap iterasi tentang proses optimasi.

    TD = fmincon(@(x)dTUCdTD(x(1), xi, D, cp, cb, ch, co, theta, beta, alpha, k, A, delta, TP, b), x0, [], [], [], [], lb, ub, [], options);

    fprintf('TD = %.4f\n', TD);
    fprintf('xi = %.4f\n', xi);
end

function F = dTUCdTD(TD, xi, D, cp, cb, ch, co, theta, beta, alpha, k, A, delta, TP, b)
    F = ((D*(ch + xi)*(delta*xi + 1)*(exp((TD*theta)/(delta*xi + 1)) - 1))/theta + ...
         (D*cp*exp(TP*alpha + TD*theta)*(b*exp(-TP*k) + 1))/(b + 1) + ...
         (D*cb*exp(TP*alpha + TD*theta)*(b*exp(-TP*k) + 1)*(exp(TP*beta) - 1))/(A*beta))/TD - ...
         (co - (D*(ch + xi)*(delta*xi + 1)*(TD*theta - exp((TD*theta)/(delta*xi + 1)) + delta*xi - delta*xi*exp((TD*theta)/(delta*xi + 1)) + 1))/theta^2 + ...
         (D*cp*exp(TP*alpha)*(b*exp(-TP*k) + 1)*(exp(TD*theta) - 1))/(theta*(b + 1)) + ...
         (D*cb*exp(TP*alpha)*(b*exp(-TP*k) + 1)*(exp(TP*beta) - 1)*(exp(TD*theta) - 1))/(A*beta*theta))/TD^2;
end
