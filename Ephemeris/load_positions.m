function [data, dates] = load_positions(fname)
% Loads data from ephemerides file
% Reads an ephemeris file generated from NASA JPL's HORIZONS tool.
% 
% INPUT
% fname - string
%         File name to load from
% 
% OUTPUT
% data  - N x 7 Double matrix
%         N rows of data
% dates - N x 1 cell
%         elements of this cell contain the string that represents the date
%         of the particular data point.
% 
% DEPENDENCIES
% 
% @author: Matt Marti
% @date: 2018-12-04

fid = fopen(fname,'r');

% Skip Header
nheader = 0;
while ~feof(fid)
    nheader = nheader + 1;
    str = fgetl(fid);
    if strcmp(str,'$$SOE')
        break;
    end
end

% Count number of lines
n = 0;
str = '';
while ~feof(fid) && ~strcmp(str, '$$EOE')
    str = fgetl(fid);
    n = n + 1;
end
n = n - 1;
fseek(fid, 0, 'bof');

% Skip Header
for i = 1:nheader
    fgetl(fid);
end

% Read data
i = 1;
data = zeros(n,7);
dates = cell(n,1);
while ~feof(fid) && i <= n
    str = fgetl(fid);
    if strcmp(str, '##EOE')
        data = data(1:i,:);
        break;
    end
    datacell = textscan(str, '%s', 'Delimiter', ',');
    d = str2double(datacell{1});
    data(i,:) = d([1, 3:8])';
    dates{i} = datacell{1}{2};
    i = i + 1;
end
fclose(fid);

end