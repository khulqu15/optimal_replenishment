function xi = menentukannilaixi()
    D = 100 * 10^6;
    theta = 0.2;
    ch = 0.002;
    delta = 850;
    TD = 0.2145;


    % Define the nonlinear function to solve
    f = @(xi) (-D * ((2 * ch * delta - exp((TD * theta)) / (delta * xi + 1)) + TD * theta + (4 * delta * xi) + (3 * delta^2 * xi^2) - (3 * delta^2 * xi^2 * exp((TD * theta) / (delta * xi + 1))) - ...
        (2 * ch * delta * exp((TD * theta) / (delta * xi + 1))) - (4 * delta * xi * exp((TD * theta) / (delta * xi + 1))) + (2 * ch * delta^2 * xi) + (2 * TD * delta * theta * xi) - (2 * ch * delta^2 * xi * exp((TD * theta) / (delta * xi + 1))) + ...
        (TD * ch * delta * theta) + (TD * ch * delta * theta * exp((TD * theta) / (delta * xi + 1))) + (TD * delta * theta * xi * exp((TD * theta) / (delta * xi + 1)) + 1))) / (TD * theta^2);

    % Set options for fsolve to display each iteration
    options = optimset('Display', 'off');

    xi0 = 0;
    
    % Solve using fmincon
    xi = fsolve(f, xi0, options);

    fprintf('Optimized xi = %.7f\n', xi);

    return;
end
