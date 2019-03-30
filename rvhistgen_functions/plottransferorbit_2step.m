function plottransferorbit_2step( mu, x_A_0, x_A_1, t_A, x_B_1, t_B, ...
    x_C, t_C, EorVflag, fig )
% Plots the 2 step transfer orbit from Earth to Venus
% "2 step" meaning there is a burn at 1/2*tof
% 
% @arg
% mu       - double
%            Solar Gravity constant
% x_A_0    - 6 x 1 double matrix
%            State at position A before burn
% x_A_1    - 6 x 1 double matrix
%            State at position A after burn
% t_A      - double
%            Time of position A burn
% x_B_1    - 6 x 1 double matrix
%            State at position B after burn
% t_B      - double
%            Time of positoin B at burn
% x_C      - 6 x 1 double matrix
%            State at position C
% t_C      - double
%            Time at position C
% EorVflag - 
% fig      - int
%            Figure Number
% 
% @author: Matt Marti
% @date: 2019-03-28


%% Calculate planet orbits

% Time of flight
thistplanetorbit = t_A + 24*3600*(0:1:400);

% Plot Earth Orbit
earthorbithist = rvhistgen_universal(mu, x_A_0, t_A, thistplanetorbit)';

% Plot Venus Orbit
venusorbithist = rvhistgen_universal(mu, x_C, t_A, thistplanetorbit)';


%% Calculate transfer orbit

% Transfer Orbit First Leg
thisttrans_AB = t_A:3600:t_B;
xhisttrans_AB = rvhistgen_universal(mu, x_A_1, t_A, thisttrans_AB)';

% Transfer Orbit First Leg
thisttrans_BC = t_B:3600:t_C;
xhisttrans_BC = rvhistgen_universal(mu, x_B_1, t_B, thisttrans_BC)';


%% Plot Planet Orbits

% Open figure
figure(fig)
hold off
plot(0,0,'k.')
hold on

if EorVflag

    % Plot Earth
    plot(earthorbithist(:,1), earthorbithist(:,2), 'b', 'linewidth', 1);
    plot(x_A_0(1), x_A_0(2), 'b.', 'Markersize', 20)

    % Plot Venus
    plot(venusorbithist(:,1), venusorbithist(:,2), 'color', [1 .5 0], ...
        'linewidth', 1);
    plot(x_C(1), x_C(2), '.', 'color', [1 .5 0], 'Markersize', 20)
else
    % Plot Earth
    plot(venusorbithist(:,1), venusorbithist(:,2), 'b', 'linewidth', 1);
    plot(x_C(1), x_C(2), 'b.', 'Markersize', 20)

    % Plot Venus
    plot(earthorbithist(:,1), earthorbithist(:,2), 'color', [1 .5 0], ...
        'linewidth', 1);
    plot(x_A_0(1), x_A_0(2), '.', 'color', [1 .5 0], 'Markersize', 20)
end

% Make grid in case of exception
axis(1e8*[-2, 2, -2, 2])
axis equal
grid on, grid minor


%% Plot transition orbit

% Plot First Transfer Leg
plot(xhisttrans_AB(:,1), xhisttrans_AB(:,2), ...
    'r', 'linewidth', 1.25);

% Plot Second Transfer Leg
plot(xhisttrans_BC(:,1), xhisttrans_BC(:,2), ...
    'color', [.9 0 .6], 'linewidth', 1.25);

% Resize graph
axis(1e8*[-2, 2, -2, 2])
axis equal
drawnow

end
