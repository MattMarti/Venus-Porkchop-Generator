function [xrvhist] = satpos(mu, ephem, t0, thist)
%SATPOS Returns the ECEF position of a GPS SV
% Returns the Cartesian coordinates position of an orbiting body using the
% given set of ephemerides. Uses the keplerian orbital parameters to
% calculate position
%
% INPUTS
% mu    - double
%         massive body gravitational parameter
% ephem - 6x1 double array
%         Vector of ephemerides that describe SV orbit
% t0    - double
%         Epoch time
% thist - 1xN double array
%         Time history at which to measure SV position. First value should
%         be the time of ephemeris
% 
% OUTPUTS
% xrv - Nx6 double array
%       SV ECEF coordinates [x, y, z, vx, vy, vz]
%
% @author: Matt Marti
% @date: 2018-11-26

nparam = 6;


%% Error Detection

% Get sizes
[me, ~] = size(ephem);
[mt, nt] = size(thist);

assert(nt == 1 || mt == 1, 'Input time is not a vector');

% Flip ephem and time if they are not column arrays
if (me ~= nparam)
    ephem = ephem';
end
if (nt == 1 && mt ~= 1)
    thist = thist';
end
[me, ne] = size(ephem);

% Assert that parameter lenght is good
assert(me == nparam, 'Incorrect ephemeris length');

% Make sure inputs will result in an output that makes sense
assert((length(thist) == 1)... % One time, multiple satellites
    || (size(ephem,2) == numel(thist)) ... % One time per satellite
    || (ne == 1), ... % One satellite and multiple times
    'mismatching time and ephemeris values');


%% Ephemeris Assignment

% Ephemeris Assignment
a        = ephem(1,:);
ecc      = ephem(2,:);
M0       = ephem(3,:);
omega    = ephem(4,:);
i0       = ephem(5,:);
OMEGA0   = ephem(6,:);


%% Position

% Calculate mean anomaly
sqrta = sqrt(a);
oneoversqrta = 1./sqrta;
sqrtmu = sqrt(mu);
n = sqrtmu.*oneoversqrta.^3;
M = M0 + n .* thist;

% Compute eccentric anomaly from mean anomaly
E = M;
facabsdelFlowlim = 1000*eps;
absdelFlowlim = facabsdelFlowlim*2*pi;
k = 0;
while k < 100
    dEdM = 1 - ecc .* cos(E);
    deltaE = (E - ecc .* sin(E) - M) ./ dEdM;
    E = E - deltaE;
    if abs(deltaE) <= max([absdelFlowlim, (facabsdelFlowlim.*abs(E))])
        break;
    end
    k = k + 1;
end
sinE = sin(E);
cosE = cos(E);

% Calculate true anomoly from eccentric anomoly
sqrtonemeccsq = sqrt(1 - ecc.^2);
sinnuonemecccosE = sqrtonemeccsq.*sinE;
cosnuonemecccosE = cosE - ecc;
nu = atan2(sinnuonemecccosE, cosnuonemecccosE);

% Argument of latitude
phi = nu + omega;

% Calculate longitude of ascending node
Omega = OMEGA0;
sinOMEGA = sin(Omega);
cosOMEGA = cos(Omega);

% Calculate orbital radius
r = a.*(1 - ecc.*cos(E));

% Perifocal position
cosu = cos(phi);
sinu = sin(phi);
x = r.*cosu;
y = r.*sinu;

% Calculate inclination
sini = sin(i0);
cosi = cos(i0);

% xhat derived from rotation matrices
xhat_1 = cosOMEGA;
xhat_2 = sinOMEGA;
xhat_3 = zeros(size(cosOMEGA));
xhat = [xhat_1; xhat_2; xhat_3];

% yhat derived from rotation matrices
yhat_1 = -cosi.*sinOMEGA;
yhat_2 = cosi.*cosOMEGA;
yhat_3 = sini;
yhat = [yhat_1; yhat_2; yhat_3];

% Orbit position
rvec = x.*xhat + y.*yhat;


%% Velocity

% Time derivatives of orbital motion
onemecccosE = 1 - ecc .* cos(E);
oneoveronemecccosE = 1./onemecccosE;
Edot = n.*oneoveronemecccosE;
Edotsqrtonemeccsq = Edot.*sqrtonemeccsq;
phidot = Edotsqrtonemeccsq.*oneoveronemecccosE;

% Time derivative of u and r
udot = phidot;
rdot = a.*ecc.*sin(E).*Edot;

% Time derivatives of secular terms
xdot = rdot.*cosu - r.*sinu.*udot;
ydot = rdot.*sinu + r.*cosu.*udot;

% Orbit velocity
vvec = xdot.*xhat + ydot.*yhat;

xrvhist = [rvec; vvec];

end