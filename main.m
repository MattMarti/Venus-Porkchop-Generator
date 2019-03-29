%% E2V_PorkChop
% Uses planetary ephemerides ([X,Y,Z,Vx,Vy,Vz] coordinates) to create a
% chart of delta v values as a function of departure date and arrival date.
% 
% This script calls a function which uses ephemerides generated by NASA JPL
% Solar System Dynamics webpage for Venus and Earth to compute the delta v
% required to travel from Earth to Venus.
% 
% Data generated at https://ssd.jpl.nasa.gov/horizons.cgi#top
% 
% Assumes two body mechanics for the spacecraft and Sun bodies.
% 
% @author: Matt Marti
% @date: 2018-03-29

addpath('rvhistgen_functions');
addpath('Ephemeris');
addpath('Test');
addpath('Lambert');
clear, clc


%% Parameters

% Ephemeris files parameters
fname_Venus = 'Ephemeris/Positions_until_2125/Venus.txt'; % Venus Ephemerides file
fname_Earth = 'Ephemeris/Positions_until_2125/Earth.txt'; % Earth Ephemerides file
% fname_Venus = 'Ephemeris/Positions_Truncated/Venus.csv'; % Venus Ephemerides file
% fname_Earth = 'Ephemeris/Positions_Truncated/Earth.csv'; % Earth Ephemerides file
EorVflag = 1; % True: Destination is venus.
days_per_index = 1; % Number of days per line of data

% How to do the delta-v calculations
use_mat_file = 0; % flag to load data directly from file to avoid math
% matfile_load_name = 'porkchop_Earth_trunc.mat';
% matfile_load_name = 'porkchop_Venus_trunc.mat';
plot_progress_flag = 0; % Whether to plot orbits during computations

% What to plot
min_days_from_2020 = 0*365; % Start date for plot
max_days_from_2020 = 100*365 + min_days_from_2020; % Number of days on plot
min_days_to_travel = 1; % Minimum travel time
max_days_to_travel = 300; % Maximum travel time


%% Calculations

constants;

% Earth to Venus or Venus to Earth Status
fprintf('Earth to Venus flag: %d\n', EorVflag);

% Data file
if ~use_mat_file
    
    % Load data files
    [data_Earth, dates_Earth] = load_positions(fname_Earth);
    [data_Venus, dates_Venus] = load_positions(fname_Venus);
    
    % Convert days in file to seconds
    data_Earth(:,1) = data_Earth(:,1) * 3600*24; % Convert to seconds
    data_Venus(:,1) = data_Venus(:,1) * 3600*24; % Convert to seconds
    
    % Run Delta-V calculations
    tic
    deltav = porkchop(data_Earth, data_Venus, plot_progress_flag, EorVflag);
    toc
    
    % Unit convertions to stuff that the save file has
    data_Earth(:,1) = data_Earth(:,1) / 3600/24; % Convert back to days
    data_Venus(:,1) = data_Venus(:,1) / 3600/24; % Convert back to days
    
    % Save results to file
    save porkchop_lastdata.mat data_Earth data_Venus deltav dates_Earth dates_Venus
else
    
    % Load results from last computation
    load(matfile_load_name);
end

% Select indeces
min_indeces_from_2020 = round(min_days_from_2020/days_per_index);
max_indeces_from_2020 = round(max_days_from_2020/days_per_index);
min_indeces_to_travel = round(min_days_to_travel/days_per_index);
max_indeces_to_travel = round(max_days_to_travel/days_per_index);


%% Output

% Obtain Earliest and Lastest date strings
if EorVflag
    depart_early = dates_Earth{min_indeces_from_2020+1};
    depart_late  = dates_Earth{max_indeces_from_2020+1};
    arrive_early = dates_Venus{min_indeces_from_2020+min_indeces_to_travel+1};
    arrive_late  = dates_Venus{max_indeces_from_2020+max_indeces_to_travel+1};
else
    depart_early = dates_Venus{min_indeces_from_2020+min_indeces_to_travel+1};
    depart_late  = dates_Venus{max_indeces_from_2020+max_indeces_to_travel+1};
    arrive_early = dates_Earth{min_indeces_from_2020+1};
    arrive_late  = dates_Earth{max_indeces_from_2020+1};
end

% Print start and stop times on plot
fprintf('Earliest Departure Date: %s\n', depart_early);
fprintf('Latest Departure Date:   %s\n', depart_late);
fprintf('Earliest Arrival Date:   %s\n', arrive_early);
fprintf('Latest Arrival Date:     %s\n', arrive_late);

% Obtain departure and arrival delta_v for selected data
x = min_indeces_from_2020:max_indeces_from_2020;
y = min_indeces_to_travel:max_indeces_to_travel;
deltav_trunc = zeros(length(x), length(y));
for i = 1:length(x)
    for j = 1:length(y)
        ii = min_indeces_from_2020+i;
        jj = ii + y(j);
        deltav_trunc(i,j) = deltav( ii, jj );
    end
end

clear deltav
if EorVflag
    save porkchop_2Venus
else
    save porkchop_2Earth
end
deltav_trunc = min(deltav_trunc, DELTAV_MAX_CONSIDER);

% Plot Delta V map
figure(1)
hold off
contourf(x,y,min(deltav_trunc, DELTAV_MAX_CONSIDER)')
colormap jet; % Set the colors of the contour
colorbar vert; % Make vertical color bar legend
xlabel(['Days since ' depart_early]);
ylabel('Travel time (days)')
title('Delta V Contour Plot')
grid on, grid minor

% % Plot Delta V map
% figure(1)
% hold off
% idum = 20*365+6:25*365;
% contourf(0:length(idum)-1,y,deltav_trunc(idum,:)')
% colormap jet; % Set the colors of the contour
% colorbar vert; % Make vertical color bar legend
% xlabel(['Days since ' dates_from{idum(1)}]);
% ylabel('Travel time (days)')
% title('Delta V Contour Plot')
% grid on, grid minor