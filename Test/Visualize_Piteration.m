%% Visualize Piteration
% 
% Script to visualize the transfer orbit after a p-iteration. Used to debug
% that the transfer orbit calculation works correctly.
% 
% Taking position data from Earth and Venus assuming a 100 day transfer
% orbit. The indeces are:
%   data_from(2500,:)
%   data_to(2600,:)


clear, clc
constants
mu_sun = 1.32712440018e11; % km^3/s^2

% Earth state
t1 = 2461348.5 * 24*3600;
x1 = [ 111483852.399111;
       97002439.2587264;
       6881.45298589766;
      -20.0934823319526;
       22.3143979982634;
      -0.000422216757128169];

% Venus State
t2 = 2461448.5 * 24*3600;
x2 = [-94044538.961384;
      -53828929.9677622;
       4703786.20010104;
       16.9999115964573;
      -30.6429909048999;
      -1.40178567990664];

% Time of flight
tof = t2 - t1;
thisttrans = (t1:3600:t2)';
thistplanetorbit = t1 + 24*3600*(0:1:365)';

% Plot Earth Orbit
earthorbithist = rvhistgen(thistplanetorbit,t1,x1,mu_sun);

% Plot Venus Orbit
venusorbithist = rvhistgen(thistplanetorbit,t1,x2,mu_sun);


%% Compute transition orbit

% Compute p-iteration
r1vec = x1(1:3);
r2vec = x2(1:3);
[v1vec, v2vec] = piteration(mu_sun, r1vec, r2vec, tof, 0, 0, 0);
v1vec(2) = v1vec(2);
x1trans = [r1vec; v1vec];
x2trans = [r2vec; v2vec];

% Orbit time history
xhisttrans = rvhistgen(thisttrans,t1,x1trans,mu_sun);


%% Display

% Open figure
figure(1)
hold off
plot(0,0,'k.')
hold on

% Plot Earth
plot(earthorbithist(:,1), earthorbithist(:,2), 'b', 'linewidth', 1);
plot(x1(1), x1(2), 'b.', 'Markersize', 20)

% Plot Venus
plot(venusorbithist(:,1), venusorbithist(:,2), 'y', 'linewidth', 1);
plot(x2(1), x2(2), 'y.', 'Markersize', 20)

% Plot Transfer
plot(xhisttrans(:,1), xhisttrans(:,2), 'r', 'linewidth', 1.25);
grid on, grid minor
axis(1e8*[-2, 2, -2, 2])


