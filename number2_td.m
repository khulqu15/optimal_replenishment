function number2_td()
    % Definisi parameter masalah
    D = 100 * 10^6;
    cp = 0.005;
    lb = [0.222, 0.0015];
    ub = [1, 1];
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

    x0 = [0.00, 0.00];

    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10);

    [x, ~] = fmincon(@(x)dTUCdTD(x(1), x(2), D, cp, cb, ch, co, theta, beta, alpha, k, A, delta, TP, b), x0, [], [], [], [], lb, ub, [], options);

    TD = x(1);
    xi = x(2);

    fprintf('TD = %.4f\n', TD);
    fprintf('xi = %.5f\n', xi);
end

function F = dTUCdTD(TD, xi, D, cp, cb, ch, co, theta, beta, alpha, k, A, delta, TP, b)
    F = ((D*(ch + xi)*(delta*xi + 1)*(exp((TD*theta)/(delta*xi + 1)) - 1))/theta + ...
         (D*cp*exp(TP*alpha + TD*theta)*(b*exp(-TP*k) + 1))/(b + 1) + ...
         (D*cb*exp(TP*alpha + TD*theta)*(b*exp(-TP*k) + 1)*(exp(TP*beta) - 1))/(A*beta))/TD - ...
         (co - (D*(ch + xi)*(delta*xi + 1)*(TD*theta - exp((TD*theta)/(delta*xi + 1)) + delta*xi - delta*xi*exp((TD*theta)/(delta*xi + 1)) + 1))/theta^2 + ...
         (D*cp*exp(TP*alpha)*(b*exp(-TP*k) + 1)*(exp(TD*theta) - 1))/(theta*(b + 1)) + ...
         (D*cb*exp(TP*alpha)*(b*exp(-TP*k) + 1)*(exp(TP*beta) - 1)*(exp(TD*theta) - 1))/(A*beta*theta))/TD^2;
end
