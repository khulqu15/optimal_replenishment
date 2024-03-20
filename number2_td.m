function [TD, xi] = optimizeTDxi()
    D = 100 * 10^6;
    cp = 0.005;
    cb = 0.0315;
    ch = 0.002;
    co = 5000;
    theta = 0.2;
    beta = 81;
    b = 69.4;
    alpha = 3;
    k = 0.12 * 365;
    A = 3200;
    delta = 850;
    TP = 0.0719;

    f = @(x) dTUCdTD(x(1), x(2), D, cp, cb, ch, co, theta, beta, alpha, k, A, delta, TP, b);

    options = optimset('Display', 'iter', 'TolFun', 1e-8, 'TolX', 1e-8);
    
    x0 = [0.0001, 0.014];
    
    [sol, ~, exitflag, output] = fsolve(f, x0, options);
    
    TD = sol(1);
    xi = sol(2);
    
    fprintf("TD = %.4f, xi = %.4f", TD, xi);
end

function F = dTUCdTD(TD, xi, D, cp, cb, ch, co, theta, beta, alpha, k, A, delta, TP, b)
    F = ((D*(ch + xi)*(delta*xi + 1)*(exp((TD*theta)/(delta*xi + 1)) - 1))/theta + ...
         (D*cp*exp(TP*alpha + TD*theta)*(b*exp(-TP*k) + 1))/(b + 1) + ...
         (D*cb*exp(TP*alpha + TD*theta)*(b*exp(-TP*k) + 1)*(exp(TP*beta) - 1))/(A*beta))/TD - ...
         (co - (D*(ch + xi)*(delta*xi + 1)*(TD*theta - exp((TD*theta)/(delta*xi + 1)) + delta*xi - delta*xi*exp((TD*theta)/(delta*xi + 1)) + 1))/theta^2 + ...
         (D*cp*exp(TP*alpha)*(b*exp(-TP*k) + 1)*(exp(TD*theta) - 1))/(theta*(b + 1)) + ...
         (D*cb*exp(TP*alpha)*(b*exp(-TP*k) + 1)*(exp(TP*beta) - 1)*(exp(TD*theta) - 1))/(A*beta*theta))/TD^2;
end
