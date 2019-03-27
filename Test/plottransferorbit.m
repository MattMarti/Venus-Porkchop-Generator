function plottransferorbit( mu, x1, t1, x2, t2, fig, ccwflag, hyperflag )
% Plots the transfer orbit from Earth to Venus
% 
% @author: Matt Marti
% @date: 2018-11-26

constants
plotflag = 0;

%% Calculate planet orbits

% Time of flight
tof = t2 - t1;
thisttrans = (t1:3600:t2)';
thistplanetorbit = t1 + 24*3600*(0:1:400)';

% Plot Earth Orbit
earthorbithist = rvhistgen(thistplanetorbit,t1,x1,mu);

% Plot Venus Orbit
venusorbithist = rvhistgen(thistplanetorbit,t1,x2,mu);


%% Calculate transfer orbit
% Compute p-iteration
r1vec = x1(1:3);
r2vec = x2(1:3);
v1vec = x1(4:6);
v2vec = x2(4:6);
[v1vec, v2vec] = piteration(mu, r1vec, r2vec, tof, ccwflag, ...
    hyperflag, plotflag);
v1vec(2) = v1vec(2);
x1trans = [r1vec; v1vec];
x2trans = [r2vec; v2vec];

% Orbit time history
xhisttrans = rvhistgen(thisttrans,t1,x1trans,mu);


%% Plot Planet Orbits

% Open figure
figure(fig)
hold off
plot(0,0,'k.')
hold on

% Plot Earth
plot(earthorbithist(:,1), earthorbithist(:,2), 'b', 'linewidth', 1);
plot(x1(1), x1(2), 'b.', 'Markersize', 20)

% Plot Venus
plot(venusorbithist(:,1), venusorbithist(:,2), 'color', [1 .5 0], ...
    'linewidth', 1);
plot(x2(1), x2(2), '.', 'color', [1 .5 0], 'Markersize', 20)

% Make grid in case of exception
grid on, grid minor
axis(1e8*[-2, 2, -2, 2])
axis equal


%% Plot transition orbit

% Plot Transfer
figure(fig)
plot(xhisttrans(:,1), xhisttrans(:,2), 'r', 'linewidth', 1.25);
grid on, grid minor
axis(1e8*[-2, 2, -2, 2])
axis equal
drawnow

end

