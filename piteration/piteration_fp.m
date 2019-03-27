function [v1vec, v2vec] = piteration_fp(mu, r1vec, r2vec, tof, ccwflag,...
    hyperflag, plotflag)
%PITERATION Computes velocity vectors at t1 and t2 given the position
% Given the t2 - t1 time period of the orbit, and position r1 and r2
% vectors, this function uses p-iteration to solve for the velocity vectors
% v1 and v2. These velocity vectors can be used to calculate the delta-v
% required for an orbit transfer.
% 
% This function implements the algorithm presented in C. D. Hall
%   http://www.dept.aoe.vt.edu/~cdhall/courses/aoe4134/Apiteration.pdf
% 
% Uses false position as the root finding method instead of Newton's 
% Method. This should make the algorithm more robust, although it will 
% make it slower.
% 
% INPUT
% mu          - double
%               Standard Gravitational Parameter for massive body
% r1vec       - 3x1 double matrix
%               Time 1 position vector
% r2vec       - 3x1 double matrix
%               Time 2 position vector
% deltanu     - double
%               True anomoly difference from time 1 to time 2 (rad)
% tof         - double
%               Time of flight (t2 - t1)
% ccwflag     - bool
%               True to set p-iteration calculation to compute
%               counter-clockwise orbit
% hyperflag   - bool
%               True to set p-iteration calculation to take a hyperbolic
%               orbit
% plotflag    - bool
%               True to plot the progress of the orbit determination
% 
% OUTPUT
% v1          - 3x1 double matrix
%               Time 1 velocity vector
% v2          - 3x1 double matrix
%               Time 2 velocity vector
% 
% DEPENDENCIES
%
% @author: Matt Marti
% @date: 2018-03-26

constants;

if size(r1vec, 1) == 1
    r1vec = r1vec';
end
if size(r2vec, 1) == 1
    r2vec = r2vec';
end

% Constants
r1 = norm(r1vec(1:3));
r2 = norm(r2vec(1:3));

% Delta nu calculation
alpha_geometry = dot([0;0;1], cross(r1vec,r2vec)); % DerAstrodynamics clockwise orbit
cosdelnu = dot(r1vec,r2vec) / (r1*r2);
delnu = acos(cosdelnu);
if ccwflag % If counter-clockwise
    if alpha_geometry < 0
        delnu = 2*pi - delnu;
    end
else % Clockwise transfer
    if alpha_geometry > 0
        delnu = 2*pi - delnu;
    end
end
sindelnu = sin(delnu);
tanhalfdelnu = tan(0.5*delnu);

% Parameter calculations
k = r1*r2*(1 - cosdelnu);
l = r1 + r2;
m = r1*r2*(1 + cosdelnu);

% First guess for p
p1 = k / ( l + sqrt(2*m) );
p2 = k / ( l - sqrt(2*m) );

% iterate for guesses of alpha
% alphaguessordervec = [0.5, 0.1:0.1:0.9, 1.0:0.1:1.5, 0, -0.9, -0.8];
if hyperflag
    alphaguessvec = (1:SPACING_FP_PITERATION:50)';
else
    alphaguessvec = (0:SPACING_FP_PITERATION:1);
end

% Compute bisection bounds
tdiffvec = zeros(length(alphaguessvec),1);
pvec = zeros(length(alphaguessvec),1);
for i = 1:length(tdiffvec)
    alpha = alphaguessvec(i);
    pvec(i) = (1-alpha)*p1 + alpha*p2;
    ti = ...
        delta_t_calc(mu, tof, r1vec, r2vec, r1, r2, sindelnu, cosdelnu, ...
        tanhalfdelnu, m, k, l, pvec(i), 0);
    tdiffvec(i) = tof - ti;
end
figure(4),plot(tdiffvec)

% Find bounds
if delnu < pi
    i = 1;
    while i < length(tdiffvec)
        if tdiffvec(i)*tdiffvec(i+1) < 0
            p_a = pvec(i);
            p_b = pvec(i+1);
            t_a = tdiffvec(i);
            t_b = tdiffvec(i+1);
            if std([t_a,t_b]) > 1e6 % Account for discontinuities
                i = i + 1;
                continue;
            end
            break;
        end
        i = i + 1;
    end
    if i == length(tdiffvec)
        v1vec = NaN;
        v2vec = NaN;
        return;
    end
else
    i = length(tdiffvec) - 1;
    while i > 0
        if tdiffvec(i)*tdiffvec(i+1) < 0
            p_a = pvec(i);
            p_b = pvec(i+1);
            t_a = tdiffvec(i);
            t_b = tdiffvec(i+1);
            if std([t_a,t_b]) > 1e6 % Account for discontinuities
                i = i + 1;
                continue;
            end
            break;
        end
        i = i - 1;
    end
    if i <= 0
        v1vec = NaN;
        v2vec = NaN;
        return;
    end
end

% Bisection algorithm
i = 0;
while i < IMAX_FP_PITERATION
    
    % Bisect
    p_m = p_a - t_a * (p_b - p_a) / (t_b - t_a);
    [tof_m, ~, f, g, fdot] ...
        = delta_t_calc(mu, tof, r1vec, r2vec, r1, r2, sindelnu, ...
        cosdelnu, tanhalfdelnu, m, k, l, p_m, plotflag);
    t_m = tof - tof_m;
    if (t_m*t_a < 0)
        t_b = t_m;
        p_b = p_m;
    else
        t_a = t_m;
        p_a = p_m;
    end
    
    % Exit condition
    if abs(t_m) < PRECISION_FP_PITERATION
        break;
    end
    
    % Iterate
    i = i + 1;
    continue;
    
    
    
    % P guess using Newton Raphson
    term1 = -0.5*g/pn;
    term2 = -1.5*a*(tn-g)*( k^2 + (2*m-l^2)*pn^2 )/(m*k*pn^2);
    if ~hyperflag
        term3 = sqrt(a^3/mu)*2*k*sindelE/(pn*(k-l*pn));
    else
        term3 = - sqrt((-a)^3/mu)*2*k*sinhdelF/(pn*(k-l*pn));
    end
    ddtdp = term1 + term2 + term3;
    deltap = deltatn / ddtdp;
    pnp1 = pn + deltap;

    % Compute time of flight with parameter update
    [tnp1, ~, f, g, fdot] ...
        = delta_t_calc(mu, tof, r1vec, r2vec, r1, r2, sindelnu, cosdelnu, tanhalfdelnu, m, k, l, pnp1, plotflag);
    deltatnp1 = tof - tnp1;

    % Correct for divergence
    kdex = 0; pk = pnp1; deltatk = deltatnp1; beta = 1;
    while kdex < KMAX_PITERATION ...
            && abs(deltatk) > abs(deltatn)

        % Adjust linearized change
        beta = 0.5 * beta;
        pk = pn + beta * deltap;

        % Compute time of flight and other parameters
        [tk, ~, f, g, fdot] ...
            = delta_t_calc(mu, tof, r1vec, r2vec, r1, r2, sindelnu, cosdelnu, tanhalfdelnu, m, k, l, pk, plotflag);
        if isnan(tk)
            kdex = kdex + 1;
            continue;
        end

        % Change in time
        deltatk = tof - tk;
        delthist = [delthist, deltatk];

        % Iterate
        kdex = kdex + 1;
    end

    % Iterate
    pn = pk;
    deltatn = deltatk;
    i = i + 1;
end

% Velocity calc using f-g method
gdot  = 1 - r1/p_m*(1-cosdelnu);
v1vec = (r2vec - f*r1vec(1:3))/g;
v2vec = fdot*r1vec(1:3) + gdot*v1vec;

% Orbit plot    
if plotflag
    try
        hold on
        plotorbit(mu, [r1vec(1:3); v1vec] , 0:24*3600:tof, ...
            2, [1,0,0] );
        axis equal
        drawnow
        figure(3), plot(delthist)
    catch
        % Intentionally left blank
    end
end

end

function [t, a, f, g, fdot, hyperflag, sindelE, sinhdelF] ...
    = delta_t_calc(mu, tof, r1vec, r2vec, r1, r2, sindelnu, cosdelnu, tanhalfdelnu, m, k, l, p, plotflag)
% Computes the time of flight based on given parameter p

% Compute semi-major axis and determine if hyperbolic
a = (m*k*p) / ( (2*m-l^2)*p^2 + 2*k*l*p - k^2 );
if a < 0
    hyperflag = 1;
else
    hyperflag = 0;
end

% Compute Gauss terms
f = 1 - r2*(1 - cosdelnu)/p;
fdot = sqrt(mu/p) * tanhalfdelnu * ((1 - cosdelnu)/p - 1/r1 - 1/r2);
g = r1*r2*sindelnu/sqrt(mu*p);

% Determine eccentric anomoly change
if ~hyperflag
    sindelE = -r1*r2*fdot/sqrt(mu*a);
    cosdelE = 1 - r1/a*(1-f);
    sinhdelF = [];
    try
        delE = atan2(sindelE, cosdelE);
    catch
        t = NaN;
        return;
    end
    t = g + sqrt(a^3/mu)*(delE - sin(delE));
else
    sindelE = [];
    coshdelF = 1 - r1/a*(1-f);
    delF = acosh(coshdelF);
    sinhdelF = sinh(delF);
    t = g + sqrt((-a)^3/mu)*(sinhdelF - delF);
end

% Orbit plot
if plotflag
    gdot  = 1 - r1/p*(1-cosdelnu);
    v1vec = (r2vec(1:3) - f*r1vec(1:3))/g;
    v2vec = fdot*r1vec(1:3) + gdot*v1vec; %#ok
    try
        hold on
        plotorbit(mu, [r1vec(1:3); v1vec] , 0:24*3600:tof, ...
            2, [1,0,0] );
        axis equal
        drawnow
    catch
        % Intentionally left blank
    end
end

end
