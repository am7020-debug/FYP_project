clear
clc

% Defining system features
system_capacity = 100; % System capacity in GW
avg_CF = 0.6; % Average capacity factor, 60%
red_CF = 0.1; % Reduced capacity factor, 10%
hours = 720; % Hours in a 30-day period

% Calculating expected and reduced energy generation in GWh
exp_energy = system_capacity * avg_CF * hours;
red_energy = system_capacity * red_CF * hours;

% Calculating overall shortage in GWh and converting to TWh
shortage_GWh = exp_energy - red_energy;
TWh_shortage = shortage_GWh / 1000;

disp(['Expected energy generation: ', num2str(exp_energy), ' GWh']);
disp(['Reduced energy generation: ', num2str(red_energy), ' GWh']);
disp(['Shortage: ', num2str(shortage_GWh), ' GWh']);
disp(['Shortage in TWh: ', num2str(TWh_shortage), ' TWh']);




