function optimization3()
    % Definisi parameter masalah
    A = 3200; % Kapasitas maksimal
    b = 69.4; % Konstanta integrasi yang mencerminkan pilihan waktu nol
    k = 0.12 * 365; % Konstanta pertumbuhan dikalikan dengan jumlah hari dalam setahun
    cp = 0.005; % Biaya pembelian per unit
    cb = 0.0315; % Biaya pemeliharaan per unit item selama periode pertumbuhan
    ch = 0.002; % Biaya penyimpanan per unit item per unit waktu selama periode konsumsi
    co = 5000; % Biaya pemesanan tetap per siklus
    D = 100 * 10^6; % Laju permintaan konstan (g/tahun)
    theta = 0.2; % Konstanta yang menentukan penyebaran kurva pertumbuhan
    beta = 81; % Konstanta pertumbuhan lain
    alpha = 3; % Konstanta pertumbuhan lain
    delta = 850; % Parameter yang mempengaruhi nilai xi dalam perhitungan thetaPrime

    % Nilai awal untuk optimisasi
    x0 = [0.07, 0.21, 0.05]; % Tebakan awal untuk TP, TD, dan xi

    % Batasan untuk nilai TP, TD, dan xi
    lb = [0, 0, 0]; % Batas bawah
    ub = [1, 1, 1]; % Batas atas
    
    % Konfigurasi opsi optimisasi digunakan untuk membuat atau memodifikasi 
    % optimisasi untuk algoritma yang digunakan MATLAB Optimization
    % Toolbox, Fungsi ini digunakan untuk menyesuaikan berbagai aspek dari 
    % proses optimasi, seperti algoritma yang digunakan, toleransi, kriteria 
    % berhenti, dan output yang ditampilkan.

    % Konfigurasi opsi optimasi menggunakan optimoptions
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ... % Mengatur algoritma yang digunakan oleh fmincon. 'sqp' berarti Sequential Quadratic Programming, yang efektif untuk non-linear optimization.
        'Display', 'iter'); % Mengatur tingkat tampilan output. 'iter' akan menampilkan informasi setiap iterasi tentang proses optimasi.

    % Menjalankan optimasi untuk mencari TP, TD, dan xi yang meminimalkan TUC
    [x, fval] = fmincon(@(x)calcTUC(x, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta), x0, [], [], [], [], lb, ub, [], options);
    
    % Ekstrak nilai optimal dari hasil optimasi
    TP = x(1); % Periode pemeliharaan optimal
    TD = x(2); % Periode konsumsi optimal
    xi = x(3); % Nilai xi optimal
    TUC = fval; % Total Unit Cost minimal
    
    % Menampilkan hasil optimasi
    fprintf('TP = %.4f\n', TP);
    fprintf('TD = %.4f\n', TD);
    fprintf('xi = %.4f\n', xi);
    fprintf('TUC = %.2f\n', TUC);
end

function TUC = calcTUC(x, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta)
    % Ekstraksi variabel dari input
    TP = x(1); % Periode pemeliharaan
    TD = x(2); % Periode konsumsi
    xi = x(3); % Variabel xi
    
    % Perhitungan thetaPrime berdasarkan formula yang diberikan
    thetaPrime = (1 - ((delta * xi) / (1 + delta * xi))) * theta;
    
    % Perhitungan Total Unit Cost (TUC) berdasarkan parameter yang diberikan
    Q = D / (theta * (1 + b)) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP)); % Kuantitas yang diproduksi
    PC = cp * Q; % Purchasing Cost
    HC = (((xi + ch) * D) / thetaPrime) * (1 / thetaPrime * exp(thetaPrime * TD) - 1 / thetaPrime - TD); % Holding Cost dengan penyesuaian xi dan thetaPrime
    BC = (cb * D) / (theta * beta * A) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP)) * (exp(beta * TP) - 1); % Breeding Cost
    OC = co; % Ordering Cost
    
    % Total Unit Cost dihitung dari total biaya dibagi dengan jumlah produk yang dikonsumsi
    TUC = (PC + BC + HC + OC) / TD;
end
