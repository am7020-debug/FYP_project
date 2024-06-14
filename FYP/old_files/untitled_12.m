% Define variables
Cap_e = 10000*700; % Energy capacity in kWh
DoD = 0.8; % Depth of Discharge
RT = 0.7; % Round-trip efficiency
self_dis = 0.005; % Self-discharge rate per year
Deg_t = 0.01; % Degradation rate per cycle
N_op = 3; % Number of operational cycles per year
Life = 50; % System lifetime in years
c_f = 0.055; % USD/kWh cost of gas
CBR = 0.5; % Gas consumption ratio (adjusted value)
r = 0.05; % Discount rate
c_inv = 1000; % Investment cost in USD (example value)

% Initialize cumulative values
C_fuel_total = 0;
E_out_total = 0;

% Loop through each year of the system's lifetime
for n = 1:Life
    % Calculate annual energy output (accounting for degradation)
    E_out_annual = Cap_e * DoD * RT * (1 - Deg_t)^(N_op * (n - 1));
    
    % Calculate annual fuel consumption
    F_annual = E_out_annual / CBR;
    
    % Calculate annual fuel cost
    C_fuel_annual = F_annual * c_f;
    
    % Discount annual fuel cost and energy output to present value
    C_fuel_total = C_fuel_total + C_fuel_annual / (1 + r)^n;
    E_out_total = E_out_total + E_out_annual / (1 + r)^n;
end

% Calculate levelized fuel cost
L_FUEL = C_fuel_total / E_out_total;

% Output the levelized fuel cost
disp(['Levelized Fuel Cost: ', num2str(L_FUEL), ' USD/kWh']);
display (Cap_e);