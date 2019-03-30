function deltav = transfer_deltav_2step(X_from, X_to, EorVflag, plotflag, fignum)
% Calculates delta v required for interplanetary transfer
% Computes delta-v required for a planar tranfser orbit first. Then uses
% that delta-v to compute the three dimensional transfer delta-v using a
% two step orbit, with manuever at 1/2*tof.
%
% INPUT
% X_from   - 7x1 double matrix
%            State vector of departure body: [t; x; y; z; vx; vy; vz]
% X_to     - 7x1 double matrix
%            State vector of arrival body: [t; x; y; z; vx; vy; vz]
% EorVflag - bool
%            True: Venus is target planet
%            False: Earth is target planet
% plotflag - bool
%            True to plot the minimum delta-v cost transfer orbit
% fignum   - int
%            Figure number for plot
% 
% OUTPUT
% deltav   - double
%            Delta V required to transfer
% 
% @author: Matt Marti
% @date: 2018-03-29

% Constants
global MU_SUN MAX_POSITION_ERROR MUEARTH REARTH SOIVENUS MUVENUS RVENUS SOIEARTH PARABOLIC_N_ITERATIONS_LAMBERT
mu = MU_SUN;

% Origin State (Solar coordinates)
t_A = X_from(1);
x_A_0 = X_from(2:7)';

% Target State (Solar coordinates)
t_C = X_to(1);
x_C_0 = X_to(2:7)';

% Partition out positoin vectors and compute time of flight
r_A = x_A_0(1:3);
v_A_0 = x_A_0(4:6);
r_C = x_C_0(1:3);
tof = t_C - t_A;
t_B = t_A + 0.5*tof;

% Compute "from" planet angular momentum, force destination planet to be
% coplanar
h = cross(r_A, v_A_0);
r_C_approx_1 = r_A(3)*r_C(1)/r_A(1) ...
    - h(2)*( r_A(1)*r_C(2) - r_C(1)*r_A(2) ) / ( h(3)*r_A(1) );
r_C_approx_2 = r_A(3)*r_C(2)/r_A(2) ...
    + h(1)*( r_A(1)*r_C(2) - r_C(1)*r_A(2) ) / ( h(3)*r_A(2) );
r_C_approx_3 = ( h(1)*r_A(3)*r_C(1) + h(2)*r_A(3)*r_C(2) ) ...
    / ( h(1)*r_A(1) + h(2)*r_A(2) );
if (r_C_approx_1 - r_C_approx_2) == 0
    r_C_approx = r_C_approx_1;
elseif (r_C_approx_1 - r_C_approx_3) == 0
    r_C_approx = r_C_approx_1;
elseif (r_C_approx_2 - r_C_approx_3) == 0
    r_C_approx = r_C_approx_2;
else
    r_C_approx = (r_C_approx_1 + r_C_approx_2 + r_C_approx_3) / 3;
end

% Do Lambert Counter-Clockwise only for approximate in-plane manuever
try
    r_A_approx = r_A(1:3);
    r_C_approx = [ r_C(1:2); r_C_approx ];
    ccwflag = 1;
    tol = 5e-3;
    K = PARABOLIC_N_ITERATIONS_LAMBERT;
    [v_A_1_approx, ~] = lambert(mu, r_A_approx, r_C_approx, ...
        tof, ccwflag, tol, K);
    if sum(isnan(v_A_1_approx))
        deltav = inf;
        return;
    end
catch
    deltav = inf;
    return;
end

% Propagate orbit from point A to point C to verify solution accuracy
xrv = rvhistgen_universal(mu, [r_A_approx; v_A_1_approx], t_A, t_C);
assert(max(abs(xrv(1:3) - r_C_approx)) < MAX_POSITION_ERROR, ...
    'Bad Lambert solution');

% Propagate approximate transfer orbit in 3-D to half time of flight
v_A_1 = [ v_A_1_approx(1:2); x_A_0(6) ]; % Includes Earth z-velocity
x_A_1 = [ r_A; v_A_1 ];
x_B_0 = rvhistgen_universal(mu, x_A_1, t_A, t_B);
r_B = x_B_0(1:3);
v_B_0 = x_B_0(4:6);

% Solve Lambert's problem for second leg of transfer
try
    r_B_0 = x_B_0(1:3);
    ccwflag = 1;
    tol = 5e-3;
    K = PARABOLIC_N_ITERATIONS_LAMBERT;
    [v_B_1, ~] = lambert(mu, r_B_0, r_C, ...
        t_B - t_A, ccwflag, tol, K);
    if sum(isnan(v_B_1))
        deltav = inf;
        return;
    end
catch
    deltav = inf;
    return;
end
x_B_1 = [ r_B; v_B_1 ];

% Propagate orbit from point B to point C to verify solution accuracy
xrv = rvhistgen_universal(mu, x_B_1, t_B, t_C);
assert(max(abs(xrv(1:3) - x_C_0(1:3))) < MAX_POSITION_ERROR, ...
    'Bad Lambert solution');

% Compute manuever delta-v
deltav_manuever = norm(v_B_1 - v_B_0);

% Calculate delta-v for departure and capture burns
v1inf = norm(v_A_1 - x_A_0(4:6))*1000;
try
    if EorVflag
    
        % Aerocapture delta-v at Venus to obtain 350 km circular parking orbit
        deltav_capture = 0.064;

        % Hyperbolic escape delta-v at Earth from 409 km circular parking orbit
        h_po = 409e3; % [m] Parking orbit altitude
        r1inf = SOIEARTH;
        Eescape = 0.5*v1inf^2 - MUEARTH/r1inf;
        v0escape = sqrt(2*(Eescape + MUEARTH/(h_po+REARTH)))*1e-3;
        v0circular = sqrt(MUEARTH/(h_po+REARTH))*1e-3;
        deltav_escape = abs(v0escape - v0circular);

        % Final delta-v calculation
        deltav = deltav_capture + deltav_escape;
    else
        % Aerocapture delta-v at Earth to enter atmosphere and land
        deltav_capture = 0.0;

        % Hyperbolic escape delta-v at Venus from 3509 km circular parking orbit
        h_po = 350e3; % [m] Parking orbit altitude
        r1inf = SOIVENUS;
        Eescape = 0.5*v1inf^2 - MUVENUS/r1inf;
        v0escape = sqrt(2*(Eescape + MUVENUS/(h_po+RVENUS)))*1e-3;
        v0circular = sqrt(MUVENUS/(h_po+RVENUS))*1e-3;
        deltav_escape = abs(v0escape - v0circular);

        % Final delta-v calculation
        deltav = deltav_capture + deltav_escape;
    end
    if deltav <= 0 || ~isreal(deltav)
        deltav = inf;
        return;
    end
catch
    deltav = inf;
    return;
end

% Combine delta-v to return results
deltav = deltav_manuever + deltav_capture + deltav_escape;

% Plot transfer orbit
if plotflag
    try
%         plottransferorbit( mu, x_A, t_A, x_C, t_C, v_A1, 1 );
        plottransferorbit_2step( mu, x_A_0, x_A_1, t_A, x_B_1, t_B, ...
            x_C_0, t_C, EorVflag, fignum )
    catch
        5;
    end
end

end