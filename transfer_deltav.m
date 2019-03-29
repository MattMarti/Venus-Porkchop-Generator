function deltav = transfer_deltav(X_from, X_to, plotflag, EorVflag)
%TRANSFER_DELTAV Calculates delta v given the two states
% Uses the Rendezvous model assuming no inclination change
%
% INPUT
% X_from   - 7x1 double matrix
%            State vector of departure body: [t; x; y; z; vx; vy; vz]
% X_to     - 7x1 double matrix
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
% @author: Matt Marti
% @date: 2018-03-29

% Constants
global MU_SUN MAX_POSITION_ERROR MUEARTH REARTH SOIVENUS MUVENUS RVENUS SOIEARTH
mu = MU_SUN;

% Check input
if nargin < 4
    EorVflag = 1;
end

% % Decide if earth or venus target
% if EorVflag
%     X_from = X_Earth;
%     X_to = X_Venus;
% else
%     X_from = X_Venus;
%     X_to = X_Earth;
% end

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

% Do Lambert Counter-Clockwise
passflag_ccw = 1;
try
    ccwflag = 1;
    tol = 5e-3;
    K = 32;
    [v1_ccw, v2_ccw] = lambert(mu, r1, r2, tof, ccwflag, tol, K);
    if sum(isnan(v2_ccw))
        v1_ccw = inf;
        v2_ccw = inf;
    end
catch
    v1_ccw = inf;
    v2_ccw = inf;
    passflag_ccw = 0;
end

% Delta-V Counter-Clockwise
deltav_from_ccw  = norm(x1(4:6) - v1_ccw);
deltav_to_ccw    = norm(x2(4:6) - v2_ccw);
deltav_ccw       = deltav_from_ccw + deltav_to_ccw;

% Try to compute p-iteration for Long method
pass_cw = 1;
try
    ccwflag = 0;
    tol = 5e-3;
    K = 32;
    [v1_cw, v2_cw] = lambert(mu, r1, r2, tof, ccwflag, tol, K);
    if sum(isnan(v2_cw))
        v1_cw = inf;
        v2_cw = inf;
    end
catch
    v1_cw = inf;
    v2_cw = inf;
    pass_cw = 0;
end

% Delta-V long way
deltav_from_cw = norm(x1(4:6) - v1_cw);
deltav_to_cw   = norm(x2(4:6) - v2_cw);
deltav_cw = deltav_from_cw + deltav_to_cw;

% If both long way and short way p-iteration threw exception
if ~pass_cw && ~passflag_ccw
    deltav = inf;
    return;
end

% Determine the mininum deltav requirement
if deltav_ccw <= deltav_cw
    v1          = v1_ccw;
    v2          = v2_ccw;
    deltav      = deltav_ccw;
    deltav_from = deltav_from_ccw; %#ok
    deltav_to   = deltav_to_ccw; %#ok
else
    v1          = v1_cw;
    v2          = v2_cw;
    deltav      = deltav_cw;
    deltav_from = deltav_from_cw; %#ok
    deltav_to   = deltav_to_cw; %#ok
end

% Check that the solution is correct
xrv = rvhistgen_universal(mu, [x1(1:3); v1], t1, t2);
assert(max(abs(xrv(1:3) - x2(1:3))) < MAX_POSITION_ERROR, ...
    'Bad Lambert solution');

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

% Plot transfer orbit
if plotflag
    try
        plottransferorbit( mu, x1, t1, x2, t2, v1, 1 );
    catch
        5;
    end
end

end