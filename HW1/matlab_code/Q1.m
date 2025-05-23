clc; clear; close all

U      = 250;     % [m/sec]
L      = 5;       % [m]
nu_30k = 1e-5;    % [m^2/sec]

Re_L     = U * L / nu_30k
delta    = 0.16 * 5 / Re_L^(1/7)
Re_delta = U * delta / nu_30k
eta_x    = 5 / Re_L^(3/4)
eta_y    = delta / Re_delta^(3/4)
N_x      = L / eta_x
N_y      = delta / eta_y
N        = N_x * N_y
N_t      = Re_delta^0.5
