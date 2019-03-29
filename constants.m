%% constants.m
% 
% Writes constant values to local workspace
global MAX_DAYS_TO_TRAVEL KMAX_PITERATION IMAX_PITERATION DELTAV_MAX_CONSIDER MAX_POSITION_ERROR IMAX_FP_PITERATION PRECISION_FP_PITERATION SPACING_FP_PITERATION IMAX_LAMBERT HALVINGMAX_LAMBERT PRECISION_LAMBERT PARABOLIC_TOLERANCE_LAMBERT MU_SUN GRAVCONST MASSEARTH REARTH MUEARTH SOIEARTH MASSVENUS RVENUS MUVENUS SOIVENUS

% Delta-V computing parameters
MAX_DAYS_TO_TRAVEL = 250; % Don't bother computing more than this
DELTAV_MAX_CONSIDER = 10; % Only display deltav values up to 50 km/s
MAX_POSITION_ERROR = 1e5;

% P-Iteration parameters
KMAX_PITERATION = 25; % p-iteratoin max halving steps
IMAX_PITERATION = 50; % p-iteration max number iterations

% Fixed-Point P-Iteration parameters
IMAX_FP_PITERATION = 100;
PRECISION_FP_PITERATION = 1e-6;
SPACING_FP_PITERATION = 1e-3;

% Lambert parameters
IMAX_LAMBERT = 50;
HALVINGMAX_LAMBERT = 25;
PRECISION_LAMBERT = 1e-6;
PARABOLIC_TOLERANCE_LAMBERT = 5e-3;

% Universe constants
MU_SUN = 1.32712440018e11; % km^3/s^2
GRAVCONST = 6.674e-11; % [m^3/kg/s^2] Gravitational Constant

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