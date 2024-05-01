function [TD, xi] = cok
    cp = 0.005;
    D = 100 * 10^6;
    theta = 0.2;
    b = 69.4;
    alpha = 3;
    k = 0.12 * 365;
    cb = 0.0315;
    beta = 81;
    A = 3200;
    ch = 0.002;
    co = 5000;
    delta = 850;
    TP = 0.0719;

    % Perbaiki ekspresi fungsi anonim dTUC_dTD_xi
    dTUC_dTD_xi = @(x) ((D*(ch + x(2))*(delta*x(2) + 1)*(exp((x(1)*theta)/(delta*x(2) + 1)) - 1))/theta + ...
    (D*cp*exp(TP*alpha + x(1)*theta)*(b*exp(-TP*k) + 1))/(b + 1) + ...
    (D*cb*exp(TP*alpha + x(1)*theta)*(b*exp(-TP*k) + 1)*(exp(TP*beta) - 1))/(A*beta))/x(1) - ...
    co - (D*(ch + x(2))*(delta*x(2) + 1)*(x(1)*theta - exp((x(1)*theta)/(delta*x(2) + 1)) + ...
    delta*x(2) - delta*x(2)*exp((x(1)*theta)/(delta*x(2) + 1)) + 1))/theta^2 + ...
    (D*cp*exp(TP*alpha)*(b*exp(-TP*k) + 1)*(exp(x(1)*theta) - 1))/(theta*(b + 1)) + ...
    (D*cb*exp(TP*alpha)*(b*exp(-TP*k) + 1)*(exp(TP*beta) - 1)*(exp(x(1)*theta) - 1))/(A*beta*theta)/x(1)^2;

    options = optimset('Display', 'iter');

    [sol, ~, exitflag] = fsolve(dTUC_dTD_xi, [0.2, 0.01], options);

    if exitflag > 0
        TD = sol(1);
        xi = sol(2);
        fprintf('TD Optimal = %.4f\n', TD);
        fprintf('xi Optimal = %.4f\n', xi);
    else
        disp('fsolve tidak konvergen ke solusi.');
        TD = [];
        xi = [];
    end
end
