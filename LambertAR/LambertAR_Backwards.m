%% Lambert_Backwards.m
% 
% Solves the lambert problem backwards, so that I can debug it thurouhgly
% 
% @author: Matt Marti
% @date: 2019-03-26

clear

% Consants
MUEARTH = 398600.5e9; % [m^3/s^2] Earth Gravity constant


%% Data
% Data
t1 = 212444596800;
t2 = 212455569600;
x1 = [ -25453237.0827383; ...
        146091344.286811; ...
       -2726.52790376544; ...
       -29.8633820023531; ...
       -5.16582224670029; ...
        0.00113552686025775 ];
x2 = [ -95115290.8270729; ...
       -51300108.8882472; ...
        4736254.85162797; ...
        16.7312947149885; ...
       -30.7922748926794; ...
       -1.38837841916347 ];
mu = 132712440018;

% Solution
deltav_from_short = 9.22971435972143;
deltav_to_short = 13.5898636020414;

% Compile input data
mu = mu; %#ok
r1vec = x1(1:3);
r2vec = x2(1:3);
tof = t2 - t1;
ccwflag = 0;
hyperflag = 0;
plotflag = 0;

% Run p-iteration
[v1vec, v2vec] = piteration(mu, r1vec, r2vec, tof, ccwflag, hyperflag, plotflag);

% Assertions
assert(abs(norm(v1vec - x1(4:6)) - deltav_from_short) < 1e-3, 'Bad v1');
assert(abs(norm(v2vec - x2(4:6)) - deltav_to_short) < 1e-3, 'Bad v2');

% Plot orbit data
try %#ok
    plottransferorbit( mu, x1, t1, x2, t2, 1, ccwflag, hyperflag );
end
title('Visual verification that the p-iteration was successful');


%% Stuff

% Compute norms
r1 = norm(r1vec);
r2 = norm(r2vec);

% Constants based on input
S = sqrt((r1 + r2)^3 / mu); % Eqn 26.2

% Compute geometry parameter
alpha = dot([0;0;1], cross(r1vec,r2vec)); % DerAstrodynamics clockwise orbit
cosdelnu = dot(r1vec,r2vec) / (r1*r2);
delnu = acos(cosdelnu);
DM = 1;
if ccwflag % If counter-clockwise
    if alpha < 0
        delnu = 2*pi - delnu;
        DM = -1;
    end
else % Clockwise transfer
    if alpha > 0
        delnu = 2*pi - delnu;
        DM = -1;
    end
end
A = DM*sqrt(r1*r2*(1+cosdelnu)); % Bate 5.3-27
% Up to here is right, I checked Bates for the A value

% tau value
tau = A / (r1 + r2);

% Compute lagrange parameters f, g, gdot
zvec = zeros(3,1);
Amat = [-r1vec, -v1vec, zvec; zvec, -v2vec, r2vec];
bmat = [-r2vec; r1vec];
fggdotvec = (Amat'*Amat)\(Amat')*bmat;
f = fggdotvec(1);
g = fggdotvec(2);
gdot = fggdotvec(3);

% Test that the lagrange parameters are correct
v1vec_test = (r2vec - f*r1vec)/g;
v2vec_test = (r2vec*gdot - r1vec)/g;
assert(max(abs(v1vec_test - v1vec)) < 1e-12, 'Wrong Lagrange parameters');
assert(max(abs(v2vec_test - v2vec)) < 1e-12, 'Wrong Lagrange parameters');

% Compute ktau
ktau_1 = 1 - r1*(1-f);
ktau_2 = 1 - r2*(1-gdot);
assert(abs(ktau_2 - ktau_1) < 1e-7, 'Bad ktau')
ktau = 0.5*(ktau_1 + ktau_2)

% Compute true tau
tau = g / ( (r1+r2)*sqrt(1-ktau) )

% Compute true k
k = ktau / tau
% assert(abs(k)<sqrt(2),'K is too big: |%.3f| < |k| = |%.3f|',sqrt(2),k)









[v1, v2] = lambert(mu, r1vec, r2vec, tof)





