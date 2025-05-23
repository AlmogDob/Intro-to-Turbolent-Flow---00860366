clc; clear; close all

L         = 3;             % [m]
h         = 0.3;           % [m]
g         = 9.81;          % [m/sec^2]
alpha_25C = 22.39e-6;      % [m^2/sec]
t_sec     = L^2/alpha_25C  % [sec]
t_day     = t_sec/60/60/24 % [day]

u_h       = sqrt(h*g)      % [m/sec]
t_heating = L/u_h          % [sec]