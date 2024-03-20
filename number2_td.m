function TD = number2_td
    TD0 = 0.00000000001;
    Cp = 0.005;
    D = 100 * 10^6;
    teta = 0.2;
    b = 69.4;
    alpa = 3;
    k = 0.12 * 365;
    Cb = 0.0315;
    betaa = 81;
    A = 3200;
    Ch = 0.002;
    CO = 5000;

    TP = 0.0719;

    dTUC_dTD = @(TD) ((D*Ch*(exp(TD*teta) - 1))/teta + ...
        (D*Cp*exp(TP*alpa)*exp(TD*teta)*(b*exp(-TP*k) + 1))/(b + 1) + ...
        (D*Cb*exp(TP*alpa)*exp(TD*teta)*(b*exp(-TP*k) + 1)*(exp(TP*betaa) - 1))/(A*betaa))/TD - ...
        (CO - (D*Ch*(TD - exp(TD*teta)/teta + 1/teta))/teta + ...
        (D*Cp*exp(TP*alpa)*(b*exp(-TP*k) + 1)*(exp(TD*teta) - 1))/(teta*(b + 1)) + ...
        (D*Cb*exp(TP*alpa)*(b*exp(-TP*k) + 1)*(exp(TP*betaa) - 1)*(exp(TD*teta) - 1))/(A*betaa*teta))/TD^2;

    options = optimset('Display', 'iter');

    TD = fsolve(dTUC_dTD, TD0, options);
end
