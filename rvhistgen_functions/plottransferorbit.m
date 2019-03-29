function plottransferorbit( mu, x1, t1, x2, t2, v1, fig )
% Plots the transfer orbit from Earth to Venus
% 
% @author: Matt Marti
% @date: 2019-03-28

constants


%% Calculate planet orbits

% Time of flight
thisttrans = t1:3600:t2;
thistplanetorbit = t1 + 24*3600*(0:1:400);

% Plot Earth Orbit
earthorbithist = rvhistgen_universal(mu, x1, t1, thistplanetorbit)';

% Plot Venus Orbit
venusorbithist = rvhistgen_universal(mu, x2, t1, thistplanetorbit)';


%% Calculate transfer orbit

% Orbit time history
xtrans = [x1(1:3); v1];
xhisttrans = rvhistgen_universal(mu, xtrans, t1, thisttrans)';


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
axis(1e8*[-2, 2, -2, 2])
axis equal
grid on, grid minor


%% Plot transition orbit

% Plot Transfer
figure(fig)
plot(xhisttrans(:,1), xhisttrans(:,2), 'r', 'linewidth', 1.25);
axis(1e8*[-2, 2, -2, 2])
axis equal
drawnow

end
