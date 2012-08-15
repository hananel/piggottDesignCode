%WORLD CONSTANTS
g = 9.81;
ni = 14.9*10^-6;            % Ambient viscosity
rho0 = 1.2;                 % [kg/m^3] reference density - used in power curve prediction
Rg = 8.31447;    % [J/(mol K)]
T0 = 288.15;    % [K]
g = 9.80655;    % [m/s^2]
L = 0.0065;     % [K/m]
M = 0.0289644;  % [kg/mol]
p0 = 101325;    % [Pa]
% actual site properties
T = 25 + 273.15;% [K]
h = 900;        %[m ASL]
p = p0*(1-L*h/T0)^(g*M/(Rg*L));
rho = p*M/(Rg*T);                 % [kg/m^3] site density