kb = 1.38064852e-23;
T = 273;
L = 1;
R_ar = 8.314 / 0.039948; 
rho = 0.00012786;
P = rho * R_ar * T;
sigma = 71e-12;

Kn = kb * T / (sqrt(2) * pi * sigma^2 * P * L)