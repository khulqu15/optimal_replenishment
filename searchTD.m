syms TP TD xi A b k cp cb ch co D theta beta alpha delta
thetaPrime = (1 - ((delta * xi) / (1 + delta * xi))) * theta;

Q = D / (theta * (1 + b)) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP));
PC = (cp * D) / (theta * (1 + b)) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP));
HC = (xi + ch) * D / thetaPrime * ((1 / thetaPrime) * exp(thetaPrime * TD) - 1 / thetaPrime - TD);
BC = (cb * D) / (theta * beta * A) * (exp(theta * TD) - 1) * exp(alpha * TP) * (1 + b * exp(-k * TP)) * (exp(beta * TP) - 1);
OC = co;

TUC = (PC + BC + HC + OC) / TD;

TD_derivative = diff(TUC, xi);

disp('xi = ');
disp(simplify(TD_derivative));
