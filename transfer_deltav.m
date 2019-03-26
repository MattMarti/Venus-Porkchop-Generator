function deltav = transfer_deltav(X_Earth, X_Venus, plotflag, EorVflag)
%TRANSFER_DELTAV Calculates delta v given the two states
% Uses the Rendezvous model assuming no inclination change
%
% INPUT
% X_Earth  - 7x1 double matrix
%            State vector of departure body: [t; x; y; z; vx; vy; vz]
% X_Venus  - 7x1 double matrix
%            State vector of arrival body: [t; x; y; z; vx; vy; vz]
% plotflag - bool
%            True to plot the minimum delta-v cost transfer orbit
% EorVflag - bool
%            True: Venus is target planet
%            False: Earth is target planet
% 
% OUTPUT
% deltav   - double
%            Delta V required to transfer
% 
% DEPENDENCIES
% piteration.m v 2018-11-26
% 
% @author: Matt Marti
% @date: 2018-03-18

if nargin < 4
    EorVflag = 1;
end

% Constants
constants
mu = MU_SUN;

% Decide if earth or venus target
if EorVflag
    X_from = X_Earth;
    X_to = X_Venus;
else
    X_from = X_Venus;
    X_to = X_Earth;
end

% Origin State (Solar coordinates)
t1 = X_from(1);
x1 = X_from(2:7)';

% Target State (Solar coordinates)
t2 = X_to(1);
x2 = X_to(2:7)';

% P-iteration for solar transfer orbit
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2(1) - t1(1);

% Debugging place
if x2(1) < -106913827 && x2(2) < 9075877
     5;
end
thetaE = atan2d(x1(2), x1(1));
if thetaE < 0, thetaE = thetaE + 360; end
thetaV = atan2d(x2(2), x2(1));
if thetaV < 0, thetaV = thetaV + 360; end

% Try to compute p-iteration for short method
passflag_short = 1;
hyperflag_short = 0;
try
    [v1_short, v2_short] = piteration(mu, r1, r2, tof, 0, 0, 0);
    if prod(isnan(v2_short))
        [v1_short, v2_short] = piteration(mu, r1, r2, tof, 0, 1, 0);
    end
catch
    try
        [v1_short, v2_short] = piteration(mu, r1, r2, tof, 0, 1, 0);
        hyperflag_short = 1;
    catch
        v1_short = inf;
        v2_short = inf;
        passflag_short = 0;
    end
end

% Delta-V short way
deltav_from_short  = norm(x1(4:6) - v1_short);
deltav_to_short    = norm(x2(4:6) - v2_short);
deltav_short       = deltav_from_short + deltav_to_short;

% Try to compute p-iteration for Long method
pass_long = 1;
hyperflag_long = 0;
try
    [v1_long, v2_long] = piteration(mu, r1, r2, tof, 1, 0, 0);
    if prod(isnan(v2_long))
        [v1_long, v2_long] = piteration(mu, r1, r2, tof, 1, 1, 0);
    end
catch
    try
        [v1_long, v2_long] = piteration(mu, r1, r2, tof, 1, 1, 0);
        hyperflag_long = 1;
    catch
        v1_long = inf;
        v2_long = inf;
        pass_long = 0;
    end
end

% Delta-V long way
deltav_from_long = norm(x1(4:6) - v1_long);
deltav_to_long   = norm(x2(4:6) - v2_long);
deltav_long = deltav_from_long + deltav_to_long;

% Debugging place to pause
if abs((thetaV - thetaE) - 180) < 1
    5;
end
if abs((thetaV - thetaE) - 225) < 1
    5;
end
if abs((thetaV - thetaE) - 178) < 1
    5;
end
if abs((thetaV - thetaE) - 5) < 1
    5;
end

% If both long way and short way p-iteration threw exception
if ~pass_long && ~passflag_short
    deltav = inf;
    return;
end

% Determine the mininum deltav requirement
longwayflag = 0;
if deltav_short <= deltav_long
    v1          = v1_short;
    v2          = v2_short;
    deltav      = deltav_short;
    deltav_from = deltav_from_short; %#ok
    deltav_to   = deltav_to_short; %#ok
else
    v1          = v1_long;
    v2          = v2_long;
    deltav      = deltav_long;
    deltav_from = deltav_from_long; %#ok
    deltav_to   = deltav_to_long; %#ok
    longwayflag = 1;
end

% Calculate delta-v costs/savings from transfer and aerocapture
v1inf = norm(v1 - x1(4:6))*1000;
try
    if EorVflag
    
        % Aerocapture delta-v at Venus to obtain 350 km circular parking orbit
        deltav2 = 0.064;

        % Hyperbolic escape delta-v at Earth from 409 km circular parking orbit
        h_po = 409e3; % [m] Parking orbit altitude
        r1inf = SOIEARTH;
        Eescape = 0.5*v1inf^2 - MUEARTH/r1inf;
        v0escape = sqrt(2*(Eescape + MUEARTH/(h_po+REARTH)))*1e-3;
        v0circular = sqrt(MUEARTH/(h_po+REARTH))*1e-3;
        deltav1 = abs(v0escape - v0circular);

        % Final delta-v calculation
        deltav = deltav2 + deltav1;
    else
        % Aerocapture delta-v at Earth to enter atmosphere and land
        deltav2 = 0.0;

        % Hyperbolic escape delta-v at Venus from 3509 km circular parking orbit
        h_po = 350e3; % [m] Parking orbit altitude
        r1inf = SOIVENUS;
        Eescape = 0.5*v1inf^2 - MUVENUS/r1inf;
        v0escape = sqrt(2*(Eescape + MUVENUS/(h_po+RVENUS)))*1e-3;
        v0circular = sqrt(MUVENUS/(h_po+RVENUS))*1e-3;
        deltav1 = abs(v0escape - v0circular);

        % Final delta-v calculation
        deltav = deltav2 + deltav1;
    end
    if deltav <= 0 || ~isreal(deltav)
        deltav = inf;
    end
catch
    deltav = inf;
end

% Determine if transfer orbit is hyperbolic
% hyperflag = hyperflag_short || hyperflag_long;
hyperflag = 0;
try
    xequinoctial = rv2equinoctial([r1;v1],mu);
    e = norm(xequinoctial(2:3));
    if e + 1e-6 >= 1
        hyperflag = 1;
    end
catch
    hyperflag = 1;
    5;
end

if ~hyperflag && plotflag
    try
        plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, 0 );
        5;
    catch
        5;
    end
end

if abs(x2(1) - -95115290.8270729) < 1e-1
    5;
end

end