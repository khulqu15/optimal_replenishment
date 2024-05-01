syms A b k cp cb ch co D theta beta alpha delta TP TD xi

thetaPrime = (1 - ((delta * xi) / (1 + delta * xi))) * theta;

Q = D / (theta * (1 + b)) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP));
PC = cp * Q;
HC = (((xi + ch) * D) / thetaPrime) * (1 / thetaPrime * exp(thetaPrime * TD) - 1 / thetaPrime - TD);
BC = (cb * D) / (theta * beta * A) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP)) * (exp(beta * TP) - 1);
OC = co;

TUC = (PC + BC + HC + OC) / TD;

dTUC_dTD = diff(TUC, TD);
dTUC_dxi = diff(TUC, xi);

disp('Turunan TUC terhadap TD:\n');
pretty(dTUC_dTD);
disp('Turunan TUC terhadap xi:\n');
disp(dTUC_dxi);
