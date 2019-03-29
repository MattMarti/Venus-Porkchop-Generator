function [xrvhist] = rvhistgen_universal(mu, xrv0vec, t0, thist, tol, K)
% Generates a time history of position, velocity using Universal Variables
% This function uses Universal Variables as described in Bate
% "Fundamentals of Astrodynamics", 1971 to compute a time history of
% position and velocity. The advantage of using this function over the
% keplerian elements or equinoctial elements is that Universal Variables
% are well behaved in both elliptical and hyperbolic orbits.
% 
% Note that this function is liable to diverge, as the numerical methods
% used to solve it are not 100% guarunteed to converge. A combination of
% Newton's Method and Halley's Method are used for root finding to speed
% convergence and to decrease the odds of a failure.
% 
% Also note that there is a bug where the solution is undefined at a time
% of flight equal to zero in a hyperbolic orbit. This can be solved by
% using series approximation at epoch time.
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
global PARABOLIC_TOLERANCE_LAMBERT IMAX_LAMBERT HALVINGMAX_LAMBERT
tolz = PARABOLIC_TOLERANCE_LAMBERT;

% Input checking
if nargin < 6
    K = 16;
end
if nargin < 5
    tol = 1e-9;
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

% Compute time of flight for this iteration
sqrtoomu = sqrt(oomu);
sqrtmuthist = sqrtmu*thist;
r0dotv0osqrtmu = r0dotv0*sqrtoomu;
[T, Tp, Tpp] = tof_func(xhist, ooa, r0, sqrtmuthist, r0dotv0osqrtmu, tolz, K);

% Solve Universal Time of Flight Equation
i = 0;
while i < IMAX_LAMBERT
    
    % Newton's Method
    deltax_1 = - T ./ Tp;
    xhist_1 = xhist + deltax_1;
    
    % Halley's Method
    deltax_2 = - (T.*Tp) ./ (Tp.^2 - 0.5*T.*Tpp);
    xhist_2 = xhist + deltax_2;
    
    % Compute time of flight using both methods
    [T_1, Tp_1, Tpp_1] = tof_func(xhist_1, ooa, r0, sqrtmuthist, r0dotv0osqrtmu, tolz, K);
    [T_2, Tp_2, Tpp_2] = tof_func(xhist_2, ooa, r0, sqrtmuthist, r0dotv0osqrtmu, tolz, K);
    
    % Step size halving to prevent overshoot
    alpha = 1;
    j = 0;
    while j < HALVINGMAX_LAMBERT ...
            && (max(abs(T_2)) >= max(abs(T)) && max(abs(T_1)) >= max(abs(T)))
        
        % Half step Newton
        xhist_1 = xhist + alpha*deltax_1;
        [T_1, Tp_1, Tpp_1] = tof_func(xhist_1, ooa, r0, sqrtmuthist, r0dotv0osqrtmu, tolz, K);
        
        % Half step Halley
        xhist_2 = xhist + alpha*deltax_2;
        [T_2, Tp_2, Tpp_2] = tof_func(xhist_2, ooa, r0, sqrtmuthist, r0dotv0osqrtmu, tolz, K);
        
        % Iterate
        alpha = 0.5*alpha;
        j = j + 1;
    end
    
    % Decide between using the two methods
    if max(abs(T_2)) <= max(abs(T_1)) % Choose Halley's Method
        T = T_2;
        Tp = Tp_2;
        Tpp = Tpp_2;
        deltax = deltax_2;
        xhist = xhist_2;
    else % Choose Newton's Method
        T = T_1;
        Tp = Tp_1;
        Tpp = Tpp_1;
        deltax = deltax_1;
        xhist = xhist_1;
    end
    
    % Break condition
    if max(abs(deltax))/max(abs(xhist)) < tol || max(abs(T)) == 0
        break;
    end
    
    % Iterate
    i = i + 1;
end
xsqhist = xhist.^2;
zhist = ooa*xsqhist;
Chist = C_func(zhist, tolz, K); % Bate 4.4-10
Shist = S_func(zhist, tolz, K); % Bate 4.4-11

% Evaluate f and g
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

function [T, Tp, Tpp] = tof_func(xhist, ooa, r0, sqrtmuthist, ...
                                 r0dotv0osqrtmu, tolz, K)
% Time of flight function

% Compute z
xsqhist = xhist.^2;
zhist = ooa*xsqhist; % 4.4-7

% Variables
Chist = C_func(zhist, tolz, K); % Bate 4.4-10
Shist = S_func(zhist, tolz, K); % Bate 4.4-11

% Time of flight
%     r0dotv0osqrtmu = r0dotv0*sqrtoomu;
sqrtmuti = r0dotv0osqrtmu.*xsqhist.*Chist ...
    + (1-r0*ooa).*xsqhist.*xhist.*Shist ...
    + r0.*xhist; % 4.4-14
T = sqrtmuti - sqrtmuthist;

% Derivative of time of flight
Tp = xsqhist.*Chist ...
    + r0dotv0osqrtmu.*xhist.*(1-zhist.*Shist) ...
    + r0.*(1-zhist.*Chist); % 4.4-17

% Second Derivative of time of flight
dzdxhist = 2*ooa*xhist;
dCdzhist = dCdz_func(zhist, Chist, Shist, tolz, K);
dSdzhist = dSdz_func(zhist, Chist, Shist, tolz, K);
dCdxhist = dCdzhist.*dzdxhist;
dSdxhist = dSdzhist.*dzdxhist;
Tpp = ...
    + 2*xhist.*Chist ...
    + xsqhist.*dCdxhist ...
    + r0dotv0osqrtmu.*(1-zhist.*Shist) ...
    - r0dotv0osqrtmu.*(dzdxhist.*Shist + zhist.*dSdxhist) ...
    - r0*(dzdxhist.*Chist + zhist.*dCdxhist);

end