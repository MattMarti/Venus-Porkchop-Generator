%% Lambert_Backwards.m
% 
% Solves the lambert problem backwards, so that I can debug it thurouhgly
% 
% @author: Matt Marti
% @date: 2019-03-26

clear

% Consants
MUEARTH = 398600.5e9; % [m^3/s^2] Earth Gravity constant
mu = MUEARTH;


%% Stuff

% Given
tof = 36000; % [s] Time of flight
r1vec = 1e3*[22592.145603; -1599.915239; -19783.950506]; % [m]
r2vec = 1e3*[1922.067697; 4054.157051; -8925.727465]; % [m]

% Truth values
v1vec_t = 1e3*[2.000652697; 0.387688615; -2.666947760]; % [m/s]
v2vec_t = 1e3*[-3.79246619; -1.77707641; 6.856814395]; % [m/s]

% Norms
r1 = norm(r1vec);
r2 = norm(r2vec);

% These I had to correct but now I know they're right
c = norm(r1vec-r2vec);
m = r1+r2+c;
n = r1+r2-c;
tau = 4 * tof * sqrt(mu / m^3);

% Psi, I checked with the solution on the pdf and it's right
psi0 = acos(dot(r1vec,r2vec) / (r1*r2));
psi = psi0;

% Start to compute errors
ecvec = (r2vec-r1vec)/c;
er1vec = r1vec/r1;
er2vec = r2vec/r2;

% Compute vc and vr #1
i = 1;
b = [v1vec_t(i); v2vec_t(i)];
A = [ecvec(i), er1vec(i); ecvec(i), - er2vec(i)];
vcvr_1 = A\b;

% Compute vc and vr #2
i = 2;
b = [v1vec_t(i); v2vec_t(i)];
A = [ecvec(i), er1vec(i); ecvec(i), - er2vec(i)];
vcvr_2 = A\b;

% Compute vc and vr #3
i = 3;
b = [v1vec_t(i); v2vec_t(i)];
A = [ecvec(i), er1vec(i); ecvec(i), - er2vec(i)];
vcvr_3 = A\b;

% Average the results together
vcvr = 1/3 * (vcvr_1 + vcvr_2 + vcvr_3);

% Compute true values of x and y parameters
A = [sqrt(n), sqrt(m); -sqrt(n), sqrt(m)];
b = vcvr*sqrt(m*n/mu);
xy = A\b;
x = xy(1);
y = xy(2);

% Compute sigma based on equation 8
sigmasq = (1 - y^2) / (1 - x^2);
sigma1 = sqrt(sigmasq);
sigma2 = sqrt(4*r1*r2/m^2 * cos(0.5*psi)^2);
assert(abs(sigma1-sigma2) < 1e-6, 'Bad sigma');


%% Lambert Equation

N = 1;

phi = @(u) acot(u/sqrt(1-u^2)) - 1/(3*u)*(2+u^2)*sqrt(1-u^2);

% Terms
omxsq = 1 - x^2;
omysq = 1 - y^2;
sqrtomxsq = sqrt(omxsq);
sqrtomysq = sqrt(omysq);

% Lambert Equation A1
F = (1/sqrtomxsq^3) * (...
    - acot(x/sqrtomxsq) ...
    + acot(y/sqrtomysq) ...
    + sqrt(x*omxsq) ...
    - sqrt(y*omysq) ...
    + N*pi ...
    - tau)

% Lambert Equation 4
F = phi(x) + phi(y) + N*pi - 

% Looks like these equations for the F don't result in zero when given the
% true x and y values. This is broken! Or I'm wrong.


