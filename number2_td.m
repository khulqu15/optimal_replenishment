function [TD_opt, xi_opt] = number2_td_iterative_approach()
    % Mendeklarasikan fungsi dengan nama `number2_td_iterative_approach` yang mengembalikan TD dan xi optimal

    % Parameter
    D = 100 * 10^6; % Jumlah permintaan
    theta = 0.2; % Parameter model
    delta = 850; % Parameter model
    Ch = 0.002; % Biaya pemeliharaan per unit
    CO = 5000; % Biaya pemesanan
    Cp = 0.005; % Biaya pembelian per unit
    b = 69.4; % Parameter terkait tingkat kelangsungan hidup
    alpa = 3; % Parameter model34cxfgyu0-
    k = 0.12 * 365; % Faktor penurunan tahunan dikonversi menjadi basis harian
    Cb = 0.0315; % Biaya penyimpanan per unit
    betaa = 81; % Parameter model
    A = 3200; % Parameter model
    TP = 0.0719; % Tingkat produksi

    % Range tebakan awal untuk TD dan xi
    TD_range = linspace(0.1, 0.3, 10); % Membuat array linear dari 0.1 hingga 0.3 dengan 10 elemen
    xi_range = linspace(0.0001, 0.001, 10); % Membuat array linear dari 0.0001 hingga 0.001 dengan 10 elemen

    % Inisialisasi
    min_error = inf; % Menginisialisasi error minimum sebagai infiniti
    TD_opt = NaN; % Inisialisasi TD optimal sebagai NaN (Not a Number)
    xi_opt = NaN; % Inisialisasi xi optimal sebagai NaN

    options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');
    % Opsi untuk solver fsolve menggunakan algoritma Levenberg-Marquardt

    % Fungsi sistem untuk fsolve
    fungsiSistem = @(v) ...
                [((D*(Ch + v(2))*(delta*v(2) + 1)*(exp((v(1)*theta)/(delta*v(2) + 1)) - 1))/theta + ...
                (D*Cp*exp(TP*alpa + v(1)*theta)*(b*exp(-TP*k) + 1))/(b + 1) + ...
                (D*Cb*exp(TP*alpa + v(1)*theta)*(b*exp(-TP*k) + 1)*(exp(TP*betaa) - 1))/(A*betaa))/v(1) - ...
                (CO - (D*(Ch + v(2))*(delta*v(2) + 1)*(v(1)*theta - exp((v(1)*theta)/(delta*v(2) + 1)) + ...
                delta*v(2) - delta*v(2)*exp((v(1)*theta)/(delta*v(2) + 1)) + 1))/theta^2 + ...
                (D*Cp*exp(TP*alpa)*(b*exp(-TP*k) + 1)*(exp(v(1)*theta) - 1))/(theta*(b + 1)) + ...
                (D*Cb*exp(TP*alpa)*(b*exp(-TP*k) + 1)*(exp(TP*betaa) - 1)*(exp(v(1)*theta) - 1))/(A*betaa*theta))/v(1)^2
                ];
    % Definisi fungsi sistem yang kompleks untuk menemukan titik nol menggunakan fsolve

    % Mencoba kombinasi tebakan awal
    for TD_guess = TD_range
        for xi_guess = xi_range
            v0 = [TD_guess, xi_guess]; % Tebakan awal
            [solusi, ~, exitflag, ~] = fsolve(fungsiSistem, v0, options);
            TD = solusi(1);
            xi = solusi(2);

            % Hitung error dari target
            error = abs(TD - 0.2) + abs(xi - 0.0); % Menghitung kesalahan total dari nilai target

            % Update solusi optimal jika error lebih kecil
            if error < min_error && exitflag > 0 % Periksa exitflag untuk solusi valid
                min_error = error;
                TD_opt = TD;
                xi_opt = xi;
            end
        end
    end

    fprintf('TD Optimal: %.4f, xi Optimal: %.5f\n', TD_opt, xi_opt);
    % Menampilkan TD dan xi optimal
end
