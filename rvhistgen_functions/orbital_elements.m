function [a, e, M0, omega, i, Omega, nu] = orbital_elements(mu, x)
%ORBITAL_ELEMENTS 6 Keplerian Elements based on state x
% Calculates the six orbtial elements when given a position and velocity
% state x of a body
% 
% INPUTS
% mu    - double 
%         Standard Gravitational Parameter for massive body (ur mom)
% x     - 6x1 double matrix
%         Orbiting body state vector: [x; y; z; vx; vy; vz]
% 
% OUTPUTS
% a     - double
%         Semi-major Axis
% e     - double
%         Eccentricity
% M0    - double
%         Mean anomoly
% omega - double
%         Argument of Perigee (rad)
% i     - double
%         Inclination (rad)
% Omega - double
%         Longitude of Ascending Node (rad)
% nu    - double
%         True Anomoly (rad)
% 
% DEPENDENCIES
% 
% @author: Matt Marti
% @date: 2018-11-08

assert(size(x, 2) == 1 && size(x, 1) == 6, 'Bad state vector size');
    
% Partition state to position and velocity
rvec = x(1:3, 1);
vvec = x(4:6, 1);
r = norm(rvec);
v = norm(vvec);

% Energy and Semi-major axis
E = norm(vvec)^2/2 - mu/norm(rvec);
a = -0.5*mu/E;

% Angular momentum and Eccentricity
h = cross(rvec,vvec);
evec = 1/mu * ( (v^2-mu/r)*rvec - dot(rvec,vvec)*vvec );
e = norm(evec);

% Angular parameters
I = [1;0;0];
J = [0;1;0]; %#ok
K = [0;0;1];
nvec = cross(K,h);
n = norm(nvec);

% Inclination
i = acos(dot(K,h)/norm(h));
if i > pi
    i = i - pi;
end

% Longitude of Ascending Node
Omega = acos(dot(I,nvec)/norm(nvec));
if nvec(2) >= 0
    if Omega < 0
        Omega = Omega + pi;
    elseif Omega > pi
        Omega = Omega - pi;
    end
elseif nvec(2) < 0
    if Omega < pi
        Omega = Omega + pi;
    elseif Omega > 360
        Omega = Omega - pi;
    end
end

% Argument of perigee
omega = acos(dot(evec,nvec)/(e*n));
if evec(3) >= 0
    if omega < 0
        omega = omega + pi;
    elseif omega > pi
        omega = omega - pi;
    end
elseif evec(3) < 0 && abs(evec(3) > 1e-5)
    if omega < pi
        omega = omega + pi;
    end
end

% True anomoly
nu = acos(dot(rvec,evec)/(norm(rvec)*e));
if dot(rvec,vvec) > 0
    if nu < 0
        nu = nu + pi;
    elseif nu > pi
        nu = nu - pi;
    end
elseif dot(rvec,vvec) < 0 
    if nu < pi
        nu = nu + pi;
    end
end

% Determine cosnu
cosnu = cos(nu);
eccpcosnu = e + cosnu;
onepecccosnu = 1 + e*cosnu;
cosE0 = eccpcosnu / onepecccosnu;

% Determine sinnu
sinnu = sin(nu);
asqrtonemesq = a*sqrt(1 - e^2);
rvec = norm(x(1:3));
sinE0 = sinnu * rvec / asqrtonemesq;

% M0
E0 = atan2(sinE0, cosE0);
M0 = E0 - e.*sinE0;
end

