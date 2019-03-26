function deltav = porkchop(data_Earth, data_Venus, plotflag, EorVflag)
%% PORKCHOP
% 
% Uses planetary ephemerides ([X,Y,Z,Vx,Vy,Vz] coordinates) to compute the
% necessary delta v to travel from a departure planet to a destination
% planet.
% 
% This script uses ephemerides generated by NASA JPL Solar System Dynamics
% webpage for planets to compute the delta v required.
% 
% Data generated at https://ssd.jpl.nasa.gov/horizons.cgi#top
% 
% Assumes two body mechanics for the spacecraft and Sun bodies.
% 
% INPUT
% data_from  - N x 7 double matrix
%              Ephemeris for departure planet
% data_to    - N x 7 double matrix
%              Ephemeris for arrival planet
% plotflag   - bool
%              whether to plot each transfer orbit as it is computed
% 
% OUTPUT
% deltav     - double matrix, nfrom x nto transfer orbit delta v matrix
% 
% DEPENDENCIES
% transfer_deltav.m v 2018-11-08
% 
% @author: Matt Marti
% @date: 2018-03-18


%% Data

% Constants
constants;

max_seconds = MAX_DAYS_TO_TRAVEL*3600*24;


%% Calculate delta v

% Preallocate data
nfrom = size(data_Earth, 1);
nto = size(data_Venus, 1);
deltav = zeros(nfrom, nto);

% Open progress bar
progbar = waitbar(0, 'Progress');

% Loop
for i = 1:nfrom
    for j = 1:nto
        if data_Earth(i,1) < data_Venus(j,1) ...
                && data_Venus(j,1) - data_Earth(i,1) <= max_seconds
            
            % Compute the delta v
            deltav(i,j) = transfer_deltav(data_Earth(i,:), data_Venus(j,:), plotflag, EorVflag);
            
            if deltav(i,j) < 10
                5;
            end
        else
            deltav(i,j) = NaN;
        end
    end
    
    % Update progress bar
    waitbar(i/nfrom,progbar);
end

% Close progress bar
close(progbar);


end