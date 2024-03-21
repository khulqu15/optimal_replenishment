function searchXI()
    D = 100 * 10^6;
    theta = 0.2;
    delta = 850;
    ch = 0.002;
    TD = 0.5;

    equation = @(xi) calc_xi(xi, D, theta, delta, ch, TD) - xi;

    xi_initial_guess = 0.015;

    options = optimoptions('fsolve', 'Display', 'iter', 'Algorithm', 'trust-region-dogleg');
    [xi_solution] = fsolve(equation, xi_initial_guess, options);

    fprintf('xi = %.4f\n', xi_solution);
end

function xi_eq = calc_xi(xi, D, theta, delta, ch, TD)
    xi_eq = -1*((D*(2*ch*delta - exp((TD*theta)/(delta*xi + 1)) + TD*theta + 4*delta*xi + ...
        3*delta^2*xi^2 - 3*delta^2*xi^2*exp((TD*theta)/(delta*xi + 1)) - ...
        2*ch*delta*exp((TD*theta)/(delta*xi + 1)) - 4*delta*xi*exp((TD*theta)/(delta*xi + 1)) + ...
        2*ch*delta^2*xi + 2*TD*delta*theta*xi - 2*ch*delta^2*xi*exp((TD*theta)/(delta*xi + 1)) + ...
        TD*ch*delta*theta + TD*ch*delta*theta*exp((TD*theta)/(delta*xi + 1)) + ...
        TD*delta*theta*xi*exp((TD*theta)/(delta*xi + 1)) + 1)) / (TD*theta^2));
end
