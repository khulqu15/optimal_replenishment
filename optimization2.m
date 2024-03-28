        function optimization3()
        % Definisi parameter masalah
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
    
        x0 = [0.7, 0.21, 0.02];
        lb = [0.0719, 0.222, 0.0015];
        ub = [1, 1, 2];
    
        options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', 'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10);
    
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
    
        thetaPrime = (1 - ((delta * xi) / (1 + delta * xi))) * theta;
        fprintf('Theta aksen: %.4f\n', thetaPrime);
        
        modTP = TP;
        modTD = TD;
    
        Q = D / (theta * (1 + b)) * (exp(theta * modTD) - 1) * exp(alpha * modTP) * (1 + b * exp(-k * modTP));
        PC = (cp * D) / (theta * (1 + b)) * (exp(theta * modTD) - 1) * exp(alpha * modTP) * (1 + b * exp(-k * modTP));
        HC = (xi + ch) * D / thetaPrime * ((1 / thetaPrime) * exp(thetaPrime * modTD) - 1 / thetaPrime - modTD);
        fprintf('HC = %.4f\n', HC);
        BC = (cb * D) / (theta * beta * A) * (exp(theta * modTD) - 1) * exp(alpha * modTP) * (1 + b * exp(-k * modTP)) * (exp(beta * modTP) - 1);
        OC = co;
    
        TUC = (PC + BC + HC + OC) / modTD;
    end

