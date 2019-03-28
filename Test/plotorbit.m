function xrvhist = plotorbit(mu, x, thist, fig, color )
% Propagates the orbit using 2 body mechanics and draws a the plot
% 
% INPUTS
% mu    - double scalar
%         Gravitational constant of massive body
% x     - 6x1 double vector
%         [x; v] state at epoch time
% thist - 1xN double vector
%         Time values of time history
% fig   - int
%         Figure number
% color - 3x1 double vector
%         RGB values between 0 to 1
% 
% OUTPUT
% xrvhist - 6xN
% 
% @author: Matt Marti
% @date: 2018-11-29

% x size checking
if size(x,2) ~= 1
    x = x';
end
if size(thist,1) ~= 1
    thist = thist';
end


%% Propagate Orbit

% Generate orbit position / velocity time history
% xrvhist = rvhistgen(thist,thist(1),x,mu)';
xrvhist = rvhistgen_sundman(mu, x, thist(1), thist)';


%% Plots

% Open figure
figure(fig)

% Plot Orbit
plot(xrvhist(1,:), xrvhist(2,:), 'color', color, ...
    'linewidth', 1.25);
hold on
plot(xrvhist(1,1), xrvhist(2,1), '.', 'color', color, 'Markersize', 20)

end

