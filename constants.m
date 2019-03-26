%% constants.m
% 
% Writes constant values to local workspace

MU_SUN = 1.32712440018e11; % km^3/s^2
MAX_DAYS_TO_TRAVEL = 400; % Don't bother computing more than this
KMAX_PITERATION = 25; % p-iteratoin max halving steps
IMAX_PITERATION = 50; % p-iteration max number iterations
DELTAV_MAX_CONSIDER = 10; % Only display deltav values up to 50 km/s


GRAVCONST = 6.674e-11; % [m^3/kg/s^2] Gravitational Constant for Earth

% Earth Constants
MASSEARTH = 5.972e24; % [kg] Earth mass
REARTH = 6378.1e3; % [m] Equitorial Radius
MUEARTH = 398600.5e9; % [m^3/s^2] Earth Gravity constant
SOIEARTH = 9.1647e8; % [m] Earth Sphere of Influence Radius

% Earth Constants
MASSVENUS = 5.972e24; % [kg] mass
RVENUS = 6051.8e3; % [m] Equitorial Radius
MUVENUS = 324860e9; % [m^3/s^2] Gravity constant
SOIVENUS = 6.1354e8; % [m] Venus Sphere of Influence Radius