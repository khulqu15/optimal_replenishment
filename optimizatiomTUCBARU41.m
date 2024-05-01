function optimizatiomTUCBARU4()
    A = 3200;
    b = 69.4;
    k = 0.12 * 365;
    cp = 0.005;
    cb = 0.0315;
    ch = 0.002;
    co = 5000;
    D = 100 * 10^6;
    theta = 0.2;
    beta = 81;
    alpha = 3;
    delta = 850;

    x0 = [0.0, 0.1, 0.0];
    lb = [0, 0, 0.000001];
    
    ub = [1, 1, 1];
    
    options = optimoptions('fmincon','Algorithm', 'sqp','Display', 'off');

    
    [x, fval] = fmincon(@(x)calcTUC(x, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta), x0, [], [], [], [], lb, ub, [], options);

    fprintf('Optimized TP = %.4f\n', x(1));
    fprintf('Optimized TD = %.4f\n', x(2));
    fprintf('Optimized xi = %.7f\n', x(3));
    fprintf('Minimum TUC = %.2f\n', fval);
end

function TUC = calcTUC(x, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta)
    TP = x(1);
    TD = x(2);
    xi = x(3);
    
    thetaPrime = (1 - ((delta * xi) / (1 + delta * xi))) * theta;
    fprintf('Theta aksen: %.4f\n', thetaPrime);
    
    Q = D / (theta * (1 + b)) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP));
    PC = cp * Q;
    HC = (((xi + ch) * D) / thetaPrime) * ((1 / thetaPrime) * exp(thetaPrime * TD) - 1 / thetaPrime - TD);
    %HC = (((xi + ch) * D) / thetaPrime*TD) * ((1 / thetaPrime) * exp(thetaPrime * TD) - 1 / thetaPrime - TD);
    BC = (cb * D) / (theta * beta * A) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP)) * (exp(beta * TP) - 1);
    OC = co;
    
    TUC = (PC + BC + HC + OC) / TD;
end