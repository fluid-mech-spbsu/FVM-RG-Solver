%% constants for Ar
gamma = 1.67;
kb = 1.38064852e-23;
T = 273;
L = 1;
R_ar = 8.314 / 0.039948; 
rho = 0.000012786;
v = 150;

%% pressure calc
P = rho * R_ar * T;

%% Re calc
mu = 20.988 * 1e-6;
nu = mu / rho;
Re = v * L / nu;

%% Mach number
a = sqrt(P * gamma / rho);
M = v/a;

%% v1
sigma = 71e-12; % WRONG !!!
Kn = kb * T / (sqrt(2) * pi * sigma^2 * P * L)

%% v2
Kn = M/Re * sqrt(gamma * pi / 2)
%% v3 ( MFP / L )
Kn = mu / P * sqrt(pi * R_ar * T) / L

