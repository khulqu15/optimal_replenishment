function optimizationTUCBARU3()
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

    x0 = [0.0, 0.0, 0.0];
    
    lb = [0, 0.2, 0.0002];
    ub = [1, 1, 1];
    
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

    [x, fval] = fmincon(@(x)calcTUC(x, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta), x0, [], [], [], [], lb, ub, [], options);

    TP = x(1);
    TD = x(2);
    xi = x(3);
    TUC = fval;

    fprintf('TP = %.4f\n', TP);
    fprintf('TD = %.4f\n', TD);
    fprintf('xi = %.4f\n', xi);
    fprintf('TUC = %.4f\n', TUC);
end

function TUC = calcTUC(x, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta)
    TP = x(1);
    TD = x(2);
    xi = x(3);

    PC = (cp * D) / (theta * (1 + b)) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP));
    thetaPrime = (theta/(1+ (delta*xi)));
    
    %fprintf('Theta aksen: %.4f\n', thetaPrime);
    
    HC = (((xi + ch) * D * (1+ (delta*xi)))/ theta * TD)*((((1+ (delta*xi))*exp((TD*theta)/((delta*xi) + 1)))/theta)-((1+ (delta*xi))/theta)-((TD*theta)/theta));
    fprintf('HC = %.4f\n', HC);
    BC = (cb * D) / (theta * beta * A) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP)) * (exp(beta * TP) - 1);
    OC = co;

    TUC = (PC + BC + HC + OC) / TD;
end