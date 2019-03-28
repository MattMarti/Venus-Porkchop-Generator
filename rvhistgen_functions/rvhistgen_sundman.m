function [xrvhist] = rvhistgen_sundman(mu, xrv0vec, t0, thist, tol, K)
% Generates a time history of position, velocity using Universal Variables
% This function uses Universal Variables as described in Bate
% "Fundamentals of Astrodynamics", 1971 to compute a time history of
% position and velocity. The advantage of using this function over the
% keplerian elements or equinoctial elements is that Universal Variables
% are well behaved in both elliptical and hyperbolic orbits.
% 
% Note that this function is liable to diverge, and I don't know why.
% 
% @arg
% mu      - double
%           2-body gravitational parameter
% xrv0vec - 6 x 1 double matrix
%           Inertial Frame position and velocity coordinates at Epoch
% t0      - double
%           Epoch time
% thist   - N x 1 double matrix
%           Delta time values measured from Epoch
% tol     - double
%           Eccentricity tolerance to use series approximation to evaluate 
%           parabolic oribts.
% K       - int
%           Number of terms in series approximation for parabolic orbits
% 
% @return
% xrvhist - 6 x N double matrix
%           Time history of position and velocity [r1; r2; r3; v1; v2; v3]
% 
% @author: Matt Marti
% @date: 2019-03-27

% Constants
constants;
tolz = PARABOLIC_TOLERANCE_LAMBERT;

% Input checking
if nargin < 6
    K = 16;
end
if nargin < 5
    tol = 1e-15;
end
assert(numel(t0) == 1, 'Incorrect size of argument ''t0''');
assert(numel(K) == 1, 'Incorrect size of argument ''K''');
assert(numel(tol) == 1, 'Incorrect size of argument ''tol''');
assert(tol > 0, 'tol must be positive');
assert(size(xrv0vec,1) == 6 && size(xrv0vec,2) == 1, ...
    'Incorrect size of argument ''r0vec''');
assert(numel(mu) == 1, 'Incorrect size of argument ''mu''');
assert(size(thist,1) == 1, ...
    'argument ''thist'' is not either single column or single row vector');

% Time of flight calculation
thist = thist - t0;

% Magnitudes of inputs
r0vec = xrv0vec(1:3);
v0vec = xrv0vec(4:6);
r0 = norm(r0vec);
v0 = norm(v0vec);

% Orbit energy
oomu = 1/mu;
E = 0.5*v0^2 - mu / r0;

% Compute inverse of semi-major axis
ooa = - 2*E/mu;

% Compute eccentricity
r0dotv0 = dot(r0vec,v0vec);
evec = ( (v0^2-mu/r0)*r0vec - r0dotv0*v0vec )*oomu;
e = norm(evec);

% Initial guess of x
sqrtmu = sqrt(mu);
if e <= 1 - tol % Elliptical orbit
    xhist = sqrtmu*thist*ooa; % Bate 4.5-10
elseif 1 + tol <= e % Hyperbolic orbit
    sqrtma = 1/sqrt(-ooa);
    numer_4511 = -2*mu*thist;
    denom_4511 = dot(r0vec,v0vec)+sign(thist)*sqrtma*sqrtmu*(1-r0*ooa)/ooa;
    xhist = sign(thist).*sqrtma.*log(numer_4511./denom_4511); % Bate 4.5-11
else % Parabolic orbit
    xhist = 0; % ?
end

% Solve Universal Time of Flight Equation
i = 0;
sqrtoomu = sqrt(oomu);
sqrtmuthist = sqrtmu*thist;
while i < IMAX_LAMBERT
    
    % Compute z
    xsqhist = xhist.^2;
    zhist = ooa*xsqhist; % 4.4-7
    
    % Variables
    Chist = C_func(zhist, tolz, K); % Bate 4.4-10
    Shist = S_func(zhist, tolz, K); % Bate 4.4-11
    
    % Time of flight
    sqrtmuti = r0dotv0.*sqrtoomu.*xsqhist.*Chist ...
        + (1-r0*ooa).*xsqhist.*xhist.*Shist ...
        + r0.*xhist; % 4.4-14
    
    % Derivative of time of flight
    sqrtmudtidx = xsqhist.*Chist ...
        + r0dotv0.*sqrtoomu.*xhist.*(1-zhist.*Shist) ...
        + r0.*(1-zhist.*Chist); % 4.4-17
    
    % Increment x
    terri = sqrtmuthist - sqrtmuti;
    deltax = sqrtmudtidx .\ terri;
    xhistnp1 = xhist + deltax;
    
    % Break condition
    if norm(deltax)/norm(xhist) < tol, break; end
    
    % Iterate
    xhist = xhistnp1;
    i = i + 1;
end
xhist = xhistnp1;

% Evaluate f and g
xsqhist = xhist.^2;
fhist = 1 - xsqhist.*Chist/r0; % Bate 4.4-31
ghist = thist - xsqhist.*xhist.*Shist*sqrtoomu; % Bate 4.4-34

% Compute r
n = length(thist);
rvechist = zeros(3,n);
for i = 1:n
    rvechist(:,i) = fhist(i)*r0vec + ghist(i)*v0vec;
end
rhist = sqrt(sum(rvechist.^2, 1));

% Evaluate fdot and gdot
fdothist = sqrtmu./(r0.*rhist).*xhist.*(zhist.*Shist - 1); % Bate 4.4-36
gdothist = 1 - xsqhist.*Chist./rhist; % Bate 4.4-35

% Compute v
vvechist = zeros(3,n);
for i = 1:n
    vvechist(:,i) = fdothist(i)*r0vec + gdothist(i)*v0vec;
end

% Assemble output
xrvhist = [rvechist; vvechist];

end

function [ Chist ] = C_func( zhist, tol, K ) % Bate 4.4-10
Chist = zeros(size(zhist));
for i = 1:length(zhist)
    z = zhist(i);
    if z <= - tol % Elliptical orbit
        Chist(i) = (1 - cos(sqrt(z))) ./ z;
    elseif tol <= z % Hyperbolic orbit
        Chist(i) = (1 - cosh(sqrt(-z))) ./ z;
    else % Parabolic orbit
        Chist(i) = 0;
        for k = 0:K
            Chist(i) = Chist(i) + (-z).^k ./ factorial(2*k+2);
        end
    end
end
end

function [ Shist ] = S_func( zhist, tol, K ) % Bate 4.4-11
Shist = zeros(size(zhist));
for i = 1:length(zhist)
    z = zhist(i);
    if z <= - tol % Elliptical orbit
        sqrtz = sqrt(z);
        Shist(i) = (sqrtz - sin(sqrtz))./(sqrtz.^3);
    elseif tol <= z % Hyperbolic orbit
        sqrtmz = sqrt(-z);
        Shist(i) = (sinh(sqrtmz) - sqrtmz)./(sqrtmz.^3);
    else % Parabolic orbit
        Shist(i) = 0;
        for k = 0:K
            Shist(i) = Shist(i) + (-z).^k ./ factorial(2*k+3);
        end
    end
end
end