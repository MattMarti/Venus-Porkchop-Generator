%% Test_piteration
%
% Computes a p-iteration for the example presented in C. D. Hall to make
% sure that the p-iteration is working correctly
% 
% This function tests the algorithm presented in C. D. Hall
%   http://www.dept.aoe.vt.edu/~cdhall/courses/aoe4134/Apiteration.pdf
% 
% Test data from 
%   http://www.aerospacengineering.net/?p=1614
% 
% @author: Matt Marti
% @date: 2018-11-08

addpath('rvhistgen_functions')
clear


%% Test 0
% Example from some website

% Given
tof = 207 * 3600*24;
r1vec = [0.473265; -0.899215; 0]; % AU
r2vec = [0.066842;  1.561256; 0.030948]; % AU
mu = 3.964016e-14; % AU^3/s^2
v1vec_true = [ 28996.2; 15232.7;  1289.2 ];
v2vec_true = [-21147.0; 3994.5;  -663.3 ];

% Call Function
[v1vec, v2vec] = piteration(mu, r1vec, r2vec, tof, 0, 0, 0);
v1vec = v1vec*149.597870e9;
v2vec = v2vec*149.597870e9;

% Assertion
p = 1e-1;
for i = 1:3
    assert(abs(v1vec_true(i) - v1vec(i)) < p, 'Bad value');
    assert(abs(v2vec_true(i) - v2vec(i)) < p, 'Bad value');
end

% Display on plot
thist_trans = 0:24*3600:tof;
figure(1), hold off, plot(0,0);
plotorbit(mu, [r1vec; v1vec/149.597870e9], thist_trans, 1, [0,0,1]);


%% Test 1
% Earth to Venus transfer orbit
% Venus slightly ahead of Earth
% Short transfer case
% Usually works

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
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2 - t1;
longwayflag = 0;
hyperflag = 0;
plotflag = 1;

% Run p-iteration
figure(2), hold off, plotorbit(mu, x1, t1:24*3600:t1+400*24*3600, 2, [0,0,1] ); plotorbit(mu, x2, t1:24*3600:t1+400*24*3600, 2, [1,0.25,0.25] );
[v1, v2] = piteration(mu, r1, r2, tof, longwayflag, hyperflag, plotflag);

% Assertions
assert(abs(norm(v1 - x1(4:6)) - deltav_from_short) < 1e-3, 'Bad v1');
assert(abs(norm(v2 - x2(4:6)) - deltav_to_short) < 1e-3, 'Bad v2');

% Plot orbit data
plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, hyperflag );

%% Test 1.2
% Earth to Venus transfer orbit
% Venus 225 degrees ahead of Earth
% Long transfer case
% Make sure doesn't break

% Data
t1 = 212444596800;
t2 = 212461876800;
x1 = [ -25453237.0827383; ...
        146091344.286811; ...
       -2726.52790376544; ...
       -29.8633820023531; ...
       -5.16582224670029; ...
        0.00113552686025775 ];
x2 = [  88225172.8216755; ...
       -61430480.2164939; ...
       -5984741.57543087; ...
        19.8810821383247; ...
        28.5235629307946; ...
       -0.756139341675039 ];
mu = 132712440018;

% Note that the correct v1 should be about [-28.4; 0.4; 0]
% Use: 
% hold off, plotorbit(mu, [x1(1:3); -27.9; -1.0; 0] , t1:24*3600:t2, 1, [1,0,0] ); plotorbit(mu, x1, t1:24*3600:t1+400*24*3600, 1, [0,0,1] ); plotorbit(mu, x2, t2:24*3600:t2+400*24*3600, 1, [1,0.5,0.5] );

% Compile input data
mu = mu; %#ok
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2 - t1;
longwayflag = 1;
hyperflag = 0;
plotflag = 1;

% Run p-iteration
clc
figure(2), hold off, plotorbit(mu, x1, t1:24*3600:t1+400*24*3600, 2, [0,0,1] ); plotorbit(mu, x2, t1:24*3600:t1+400*24*3600, 2, [1,0.25,0.25] );
[v1, v2] = piteration(mu, r1, r2, tof, longwayflag, hyperflag, plotflag);

% Plot orbit data
plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, hyperflag );

% Assertion
assert(~sum(isnan(v1)), 'bad velocity')
assert(~sum(isnan(v2)), 'bad velocity')
assert(abs(v1(1) - -27.9) < 1, 'bad velocity')
assert(abs(v1(2) - -1.0) < 1, 'bad velocity')
assert(abs(v1(3) - 0) < 1, 'bad velocity')



%% Test 2
% Earth to Venus transfer orbit
% Venus 225 degrees ahead of Earth
% Long transfer case
% Make sure doesn't break

% Data
t1 = 212444596800;
t2 = 212461876800;
x1 = [ -25453237.0827383; ...
        146091344.286811; ...
       -2726.52790376544; ...
       -29.8633820023531; ...
       -5.16582224670029; ...
        0.00113552686025775 ];
x2 = [  88225172.8216755; ...
       -61430480.2164939; ...
       -5984741.57543087; ...
        19.8810821383247; ...
        28.5235629307946; ...
       -0.756139341675039 ];
mu = 132712440018;

% Note that the correct v1 should be about [-28.4; 0.4; 0]
% Use: 
% hold off, plotorbit(mu, [x1(1:3); -27.9; -1.0; 0] , t1:24*3600:t2, 1, [1,0,0] ); plotorbit(mu, x1, t1:24*3600:t1+400*24*3600, 1, [0,0,1] ); plotorbit(mu, x2, t2:24*3600:t2+400*24*3600, 1, [1,0.5,0.5] );

% Compile input data
mu = mu; %#ok
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2 - t1;
longwayflag = 1;
hyperflag = 0;
plotflag = 1;

% Run p-iteration
clc
figure(2), hold off, plotorbit(mu, x1, t1:24*3600:t1+400*24*3600, 2, [0,0,1] ); plotorbit(mu, x2, t1:24*3600:t1+400*24*3600, 2, [1,0.25,0.25] );
[v1, v2] = piteration(mu, r1, r2, tof, longwayflag, hyperflag, plotflag);

% Plot orbit data
plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, hyperflag );

% Assertion
assert(~sum(isnan(v1)), 'bad velocity')
assert(~sum(isnan(v2)), 'bad velocity')
assert(abs(v1(1) - -27.9) < 1, 'bad velocity')
assert(abs(v1(2) - -1.0) < 1, 'bad velocity')
assert(abs(v1(3) - 0) < 1, 'bad velocity')


%% Test 3
% Earth to Venus transfer orbit
% Venus 178 degrees ahead of Earth
% Short transfer case
% Make sure doesn't break

% Data
t1 = 212444596800;
t2 = 212459284800;
x1 = [ -25453237.0827383; ...
        146091344.286811; ...
       -2726.52790376544; ...
       -29.8633820023531; ...
       -5.16582224670029; ...
        0.00113552686025775 ];
x2 = [ 13436858.421996; ...
      -106816644.171782; ...
      -2290985.11979151; ...
       34.47205747715; ...
       4.4451082242711; ...
      -1.92859237476142 ];
mu = 132712440018;

% Compile input data
mu = mu; %#ok
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2 - t1;
longwayflag = 0;
hyperflag = 0;
plotflag = 0;

% Run p-iteration
[v1, v2] = piteration(mu, r1, r2, tof, longwayflag, hyperflag, plotflag);

% Assertion
assert(~sum(isnan(v1)), 'bad velocity')
assert(~sum(isnan(v2)), 'bad velocity')

% Plot orbit data
plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, hyperflag );


%% Test 4
% Earth to Venus transfer orbit
% Venus 180 degrees ahead of Earth
% Short transfer case
% Make doesn't break

% Data
t1 = 212444596800;
t2 = 212459457600;
x1 = [ -25453237.0827383; ...
        146091344.286811; ...
       -2726.52790376544; ...
       -29.8633820023531; ...
       -5.16582224670029; ...
        0.00113552686025775 ];
x2 = [ 19368728.7896743; ...
      -105883139.173238; ...
      -2620539.12438132; ...
       34.1662808863502; ...
       6.35638262248577; ...
      -1.8847169265181 ];
mu = 132712440018;

% Compile input data
mu = mu; %#ok
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2 - t1;
longwayflag = 0;
hyperflag = 0;
plotflag = 0;

% Run p-iteration
[v1, v2] = piteration(mu, r1, r2, tof, longwayflag, hyperflag, plotflag);

% Assertion
assert(~sum(isnan(v1)), 'bad velocity')
assert(~sum(isnan(v2)), 'bad velocity')

% Plot orbit data
plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, hyperflag );


%% Test 4
% Earth to Venus transfer orbit
% Venus 180 degrees ahead of Earth
% Short transfer case
% Make sure equal to Homman Transfer if 2-D case

% Data
t1 = 212444596800;
t2 = 212459457600;
x1 = [ -25453237.0827383; ...
        146091344.286811; ...
        0; ...
       -29.8633820023531; ...
       -5.16582224670029; ...
        0 ];
x2 = [  19368728.7896743; ...
       -105883139.173238; ...
        0; ...
        34.1662808863502; ...
        6.35638262248577; ...
        0 ];
mu = 132712440018;

% Solution
deltav_from_short = 0; % TODO
deltav_to_short = 0; % TODO

% Compile input data
mu = mu; %#ok
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2 - t1;
longwayflag = 0;
hyperflag = 0;
plotflag = 0;

% Run p-iteration
[v1, v2] = piteration(mu, r1, r2, tof, longwayflag, hyperflag, plotflag);

% Assertion
assert(~sum(isnan(v1)), 'bad velocity')
assert(~sum(isnan(v2)), 'bad velocity')

% Assertions
assert(abs(norm(v1 - x1(4:6)) - deltav_from_short) < 1e-3, 'Bad v1');
assert(abs(norm(v2 - x2(4:6)) - deltav_to_short) < 1e-3, 'Bad v2');

% Plot orbit data
plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, hyperflag );


%% Test 5
% Earth to Venus transfer orbit
% Venus 225 degrees ahead of Earth
% Long transfer case
% Make sure doesn't break

% Data
t1 = 212444596800;
t2 = 212461876800;
x1 = [ -25453237.0827383; ...
        146091344.286811; ...
       -2726.52790376544; ...
       -29.8633820023531; ...
       -5.16582224670029; ...
        0.00113552686025775 ];
x2 = [  88225172.8216755; ...
       -61430480.2164939; ...
       -5984741.57543087; ...
        19.8810821383247; ...
        28.5235629307946; ...
       -0.756139341675039 ];
mu = 132712440018;

% Compile input data
mu = mu; %#ok
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2 - t1;
longwayflag = 0;
hyperflag = 0;
plotflag = 0;

% Run p-iteration
[v1, v2] = piteration(mu, r1, r2, tof, longwayflag, hyperflag, plotflag);

% Assertion
assert(~sum(isnan(v1)), 'bad velocity')
assert(~sum(isnan(v2)), 'bad velocity')

% Plot orbit data
plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, hyperflag );


%% Test 6
% Earth to Venus transfer orbit
% Venus 178 degrees ahead of Earth
% Long transfer case
% Make sure doesn't break

% Data
t1 = 212444596800;
t2 = 212459284800;
x1 = [ -25453237.0827383; ...
        146091344.286811; ...
       -2726.52790376544; ...
       -29.8633820023531; ...
       -5.16582224670029; ...
        0.00113552686025775 ];
x2 = [  13436858.421996; ...
       -106816644.171782; ...
       -2290985.11979151; ...
        34.47205747715; ...
        4.4451082242711; ...
       -1.92859237476142 ];
mu = 132712440018;

% Compile input data
mu = mu; %#ok
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2 - t1;
longwayflag = 0;
hyperflag = 0;
plotflag = 0;

% Run p-iteration
[v1, v2] = piteration(mu, r1, r2, tof, longwayflag, hyperflag, plotflag);

% Assertion
assert(~sum(isnan(v1)), 'bad velocity')
assert(~sum(isnan(v2)), 'bad velocity')

% Plot orbit data
plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, hyperflag );


%% Test 7
% Earth to Venus transfer orbit
% Venus 5 degrees ahead of Earth
% Short transfer case, hyperbolic orbit
% Make sure doesn't break

% Data
t1 = 212444596800;
t2 = 212450040000;
x1 = [ -25453237.0827383; ...
        146091344.286811; ...
       -2726.52790376544; ...
       -29.8633820023531; ...
       -5.16582224670029; ...
        0.00113552686025775 ];
x2 = [ -29009686.8928911; ...
        104790872.334466; ...
        3065270.38126867; ...
       -33.9125780176868; ...
       -9.43799203793275; ...
        1.82712547386842 ];
mu = 132712440018;

% Compile input data
mu = mu; %#ok
r1 = x1(1:3);
r2 = x2(1:3);
tof = t2 - t1;
longwayflag = 0;
hyperflag = 0;
plotflag = 0;

% Run p-iteration
[v1, v2] = piteration(mu, r1, r2, tof, longwayflag, hyperflag, plotflag);

% Assertion
assert(~sum(isnan(v1)), 'bad velocity')
assert(~sum(isnan(v2)), 'bad velocity')

% Plot orbit data
plottransferorbit( mu, x1, t1, x2, t2, 1, longwayflag, hyperflag );


%% Output
fprintf('PASSED: Test_piteration\n');