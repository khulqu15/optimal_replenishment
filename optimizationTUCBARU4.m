function optimizatiomTUCBARU4()
    % Inisialisasi variabel yang akan digunakan dalam kalkulasi biaya total
    A = 3200; % Area yang tersedia
    b = 69.4; % Faktor b yang digunakan dalam penghitungan Q
    k = 0.12 * 365; % Tingkat pembusukan per tahun, dari per hari menjadi per tahun
    cp = 0.005; % Biaya per unit produk
    cb = 0.0315; % Biaya pemesanan per unit produk
    ch = 0.002; % Biaya penyimpanan per unit produk per tahun
    co = 5000; % Biaya operasional tetap
    D = 100 * 10^6; % Permintaan tahunan
    theta = 0.2; % Parameter theta yang digunakan dalam penghitungan Q
    beta = 81; % Parameter beta yang digunakan dalam penghitungan biaya pemesanan
    alpha = 3; % Parameter alpha yang digunakan dalam penghitungan Q
    delta = 850; % Parameter delta yang mempengaruhi perhitungan thetaPrime

    % Definisi titik awal untuk proses optimisasi
    x0 = [0.0, 0.1, 0.0]; % Nilai awal untuk TP, TD, xi
    lb = [0, 0, 0.000001]; % Batas bawah untuk variabel TP, TD, xi
    ub = [1, 1, 1]; % Batas atas untuk variabel TP, TD, xi
    
    % Konfigurasi opsi untuk algoritma optimisasi menggunakan Sequential Quadratic Programming (SQP)
    % Algoritma SQP dipilih karena kemampuannya dalam menangani masalah optimisasi non-linear dengan efisien.
    % SQP sangat cocok untuk masalah dengan banyak batasan dan variabel non-linear.
    % 'Display', 'off' dipilih untuk menghindari cetak output iterasi pada window command, mempercepat proses optimisasi.
    options = optimoptions('fmincon','Algorithm', 'sqp','Display', 'off');
    
    % Memanggil fungsi fmincon dengan fungsi tujuan calcTUC dan batasan variabel
    [x, fval] = fmincon(@(x)calcTUC(x, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta), x0, [], [], [], [], lb, ub, [], options);

    % Cetak hasil dari optimasi untuk variabel TP, TD, dan xi serta nilai minimum TUC yang didapatkan
    fprintf('Optimized TP = %.4f\n', x(1));
    fprintf('Optimized TD = %.4f\n', x(2));
    fprintf('Optimized xi = %.7f\n', x(3));
    fprintf('Minimum TUC = %.2f\n', fval);
end

function TUC = calcTUC(x, A, b, k, cp, cb, ch, co, D, theta, beta, alpha, delta)
    % Ekstraksi nilai TP, TD, dan xi dari vektor x
    TP = x(1);
    TD = x(2);
    xi = x(3);
    
    % Menghitung nilai theta' yang merupakan fungsi dari xi dan delta
    thetaPrime = (1 - ((delta * xi) / (1 + delta * xi))) * theta;
    fprintf('Theta aksen: %.4f\n', thetaPrime);
    
    % Menghitung kuantitas produk yang dipesan dengan mempertimbangkan faktor-faktor yang terlibat
    Q = D / (theta * (1 + b)) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP));
    PC = cp * Q; % Perhitungan biaya produk berdasarkan kuantitas dan biaya per unit
    
    % Perhitungan biaya penyimpanan berdasarkan kuantitas produk, biaya per unit, dan durasi penyimpanan
    HC = (((xi + ch) * D) / (thetaPrime * TD)) * ((1 / thetaPrime) * exp(thetaPrime * TD) - 1 / thetaPrime - TD);
    
    % Perhitungan biaya pemesanan berdasarkan kuantitas produk, biaya pemesanan per unit, dan faktor-faktor lainnya
    BC = (cb * D) / (theta * beta * A) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP)) * (exp(beta * TP) - 1);
    
    % Menambahkan biaya operasional tetap
    OC = co;
    
    % Menghitung Total Unit Cost (TUC) dengan membagi total biaya dengan durasi penyimpanan
    TUC = (PC + BC + HC + OC) / TD;
end
