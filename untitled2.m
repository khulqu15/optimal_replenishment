function xiOptimization()
    % Inisialisasi parameter yang dibutuhkan
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
    TP = 0.0719;
    TD = 0.6179; % Nilai TD yang diberikan

    % Nilai awal untuk xi
    xi_initial = 0.000001;

    % Menggunakan fsolve untuk mencari xi yang mengoptimalkan turunan dari TUC
    options = optimset('Display', 'iter', 'TolFun', 1e-12, 'TolX', 1e-12, 'MaxFunEvals', 10000, 'MaxIter', 10000);
    xi_solution = fsolve(@(xi) derivativeTUC(xi, TP, TD, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta), xi_initial, options);

    % Menampilkan solusi xi yang dihasilkan
    fprintf('Optimized xi = %.10f\n', xi_solution);
end

function dTUC_dxi = derivativeTUC(xi, TP, TD, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta)
    % Menghitung theta' yang merupakan fungsi dari xi dan delta
    theta = (1 - ((delta * xi) / (1 + delta * xi))) * theta;

    % Menghitung turunan HC terhadap xi
    derivativeHC = abs(-D * (2 * ch * delta - exp((TD * theta) / (delta * xi + 1)) + TD * theta + 4 * delta * xi + 3 * delta^2 * xi^2 - 3 * delta^2 * xi^2 * exp((TD * theta) / (delta * xi + 1)) - ...
        2 * ch * delta * exp((TD * theta) / (delta * xi + 1)) - 4 * delta * xi * exp((TD * theta) / (delta * xi + 1)) + 2 * ch * delta^2 * xi + 2 * TD * delta * theta * xi - 2 * ch * delta^2 * xi * exp((TD * theta) / (delta * xi + 1)) + ...
        TD * ch * delta * theta + TD * ch * delta * theta * exp((TD * theta) / (delta * xi + 1)) + TD * delta * theta * xi * exp((TD * theta) / (delta * xi + 1)) + 1) / (TD * theta^2));

    % Mengabaikan turunan PC, BC, dan OC karena tidak tergantung pada xi
    % Total turunan dari TUC
    dTUC_dxi = derivativeHC;
end
