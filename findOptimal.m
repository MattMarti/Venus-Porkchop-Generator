function [sinceStart, travelt, dv] ...
    = findOptimal(sinceStart, travelt, deltav)
% Computes the minimum delta-V transfer times and delta-v
% Uses an initial estimate obtained from the porkchop plot to search for
% the minimum delta-v of the transfer orbit.
% 
% @args
% sinceStart   - int
%                Days since start of Delta-V mesh
% travelt      - int
%                Travel time
% deltav       - double matrix
%                Precomputed Delta-V mesh
% 
% @return
% sinceStart   - int
%                Days since start of Delta-V mesh
% travelt      - int
%                Travel time
% dv           - double
%                Minimum delta-v
% 
% @author: Matt Marti
% @date: 2019-02-04

% Loop parameters
maxiter = 1000;

% Starting indeces
i = sinceStart;
j = travelt;

% Loop
k = 1;
while 1
    
    % Max iterations
    if k > maxiter
        break;
    end
    k = k + 1;
    
    % Look horizontal
    if i - 1 < 0
        % Intentionally left blank
    elseif i + 1 > size(deltav,1)
        % Intentionally left blank
    elseif deltav(i-1,j) <= deltav(i,j)
        i = i - 1;
    elseif deltav(i+1,j) <= deltav(i,j)
        i = i + 1;
    end
    
    % Look horizontal
    if j - 1 < 0
        % Intentionally left blank
    elseif j + 1 > size(deltav,2)
        % Intentionally left blank
    elseif deltav(i,j-1) <= deltav(i,j)
        j = j - 1;
    elseif deltav(i,j+1) <= deltav(i,j)
        j = j + 1;
    end
end

% Assign output
sinceStart = i;
travelt = j;
dv = deltav(i,j);

% Print output
fprintf('Days Since:  %.0f\n', sinceStart);
fprintf('Travel Time: %.0f\n', travelt);
fprintf('Delta-V:     %.3f\n', dv);

end

