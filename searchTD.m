syms A b k cp cb ch co D theta beta alpha TD TP

Q = D / (theta * (1 + b)) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP));
PC = (cp * D) / (theta * (1 + b)) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP));
BC = (cb * D) / (theta * beta * A) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP)) * (exp(beta * TP) - 1);
HC = ch * (D / theta) * ((1 / theta) * exp(theta * TD) - (1 / theta) - TD);
OC = co;

TUC = (PC + BC + HC + OC) / TD;

dTUC_dTD = diff(TUC, TD);

disp('dTUC/dTD = ');
disp(simplify(dTUC_dTD));
