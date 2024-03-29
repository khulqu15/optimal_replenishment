function searchXI()
    lb = 0.0015;
    ub = 1;
    D = 100 * 10^6;
    theta = 0.2;
    delta = 850;
    ch = 0.002;
    TD = 0.222;

    objective = @(xi) calc_xi(xi, D, theta, delta, ch, TD);

    xi0 = 0.00;

    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter');

    [xi_solution, fval] = fmincon(objective, xi0, [], [], [], [], lb, ub, [], options);

    fprintf('xi = %.4f\n', xi_solution);
end

function xi_eq = calc_xi(xi, D, theta, delta, ch, TD)
    xi_eq = ((D*(2*ch*delta - exp((TD*theta)/(delta*xi + 1)) + TD*theta + 4*delta*xi + ...
        3*delta^2*xi^2 - 3*delta^2*xi^2*exp((TD*theta)/(delta*xi + 1)) - ...
        2*ch*delta*exp((TD*theta)/(delta*xi + 1)) - 4*delta*xi*exp((TD*theta)/(delta*xi + 1)) + ...
        2*ch*delta^2*xi + 2*TD*delta*theta*xi - 2*ch*delta^2*xi*exp((TD*theta)/(delta*xi + 1)) + ...
        TD*ch*delta*theta + TD*ch*delta*theta*exp((TD*theta)/(delta*xi + 1)) + ...
        TD*delta*theta*xi*exp((TD*theta)/(delta*xi + 1)) + 1)) / (TD*theta^2))^2;
end
