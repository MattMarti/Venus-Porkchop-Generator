%% Test_rvhistgen_universal
% Test the rvhistgen function that uses universal variables
% 
% @author: Matt Marti
% @date: 2019-03-27

clear


%% Test 1

% Inputs
mu = 1;
t = 2;
xrv0vec = [1; 0; 0; 0; 0; 1.1];

% True Outputs
xrv_true = [-0.321; 0; 1.236; -0.8801; 0; -0.039 ];

% Run function
[xrv] = rvhistgen_universal(mu, xrv0vec, 0, t);

% Compare outputs
diff = xrv - xrv_true;
for i = 1:6
    assert(abs(diff(i)) < 2e-3, 'Bad value');
end


%% Test 2
% Multiple times

% Inputs
mu = 1;
xrv0vec = [1; 0; 0; 0; 0; 1.1];

% Time of flight
t1 = 0;
t2 = 18;
thist = (t1:.1:t2);

% Plot Earth Orbit
[a, e, M0, omega, i, Omega, ~] = orbital_elements(mu, xrv0vec);
ephem = [a, e, M0, omega, i, Omega];
xrvtrue = satpos(mu, ephem, t1, thist);

% Run function
[xrv] = rvhistgen_universal(mu, xrv0vec, t1, thist);

% Compare outputs
diff = xrv - xrvtrue;
assert(max(max(abs(diff))) < 1e-10, 'Bad value');


% %% Test 3
% % Time history of values
% 
% constants
% mu = 1.32712440018e11; % km^3/s^2
% 
% % Earth state
% t1 = 2461348.5 * 24*3600;
% x1 = [ 111483852.399111;
%        97002439.2587264;
%        6881.45298589766;
%       -20.0934823319526;
%        22.3143979982634;
%       -0.000422216757128169];
% t2 = t1 + 180*24*3600;
% 
% % Plot Earth Orbit
% [a, e, M0, omega, i, Omega, ~] = orbital_elements(mu, x1);
% ephem = [a, e, M0, omega, i, Omega];
% xrvtrue = satpos(mu, ephem, t1, t2);
% 
% % Run function
% [xrv] = rvhistgen_universal(mu, x1, t1, t2);
% 
% % Compare outputs
% diff = xrv - xrvtrue;
% assert(max(max(abs(diff))) < 1e-3, 'Bad value');
% 
% 
% %% Test 4
% % Time history of values
% 
% constants
% mu = 1.32712440018e11; % km^3/s^2
% 
% % Earth state
% t1 = 2461348.5 * 24*3600;
% x1 = [ 111483852.399111;
%        97002439.2587264;
%        6881.45298589766;
%       -20.0934823319526;
%        22.3143979982634;
%       -0.000422216757128169];
% t2 = t1 + 180*24*3600;
% 
% % Time of flight
% tof = t2 - t1;
% thist = (t1:3600:t2);
% 
% % Plot Earth Orbit
% [a, e, M0, omega, i, Omega, ~] = orbital_elements(mu, x1);
% ephem = [a, e, M0, omega, i, Omega];
% xrvtrue = satpos(mu, ephem, t1, t2);
% 
% % Run function
% [xrv] = rvhistgen_universal(mu, x1, 0, thist);
% 
% % Compare outputs
% diff = xrv - xrvtrue;
% assert(max(max(abs(diff))) < 5e-2, 'Bad value');


%% Output
fprintf('PASSED: Test_rvhistgen_universal\n');