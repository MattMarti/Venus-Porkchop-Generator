function [ dSdzhist ] = dSdz_func( zhist, Chist, Shist, tol, K )
% Derivative of S parameter with respect to z
% Implementation of the equation Bate 5.3-20 and 5.3-25 from "Fundamentals 
% of Astrodynamics", 1971.
% 
% @arg
% zhist    - 1 x N double matrix
%            z parameter value
% Chist    - 1 x N double matrix
%            C parameter value
% Shist    - 1 x N double matrix
%            S parameter value
% tol      - double
%            Tolerance to apply the series approximation for parabolic orbit.
% K        - int
%            Series approximation number of iterations
% 
% @return
% dSdzhist - 1 x N double matrix
%            Derivative of "Sine" parameter
% 
% @author: Matt Marti
% @date: 2019-03-28

% Preallocate
dSdzhist = zeros(size(zhist));

% Loop
for i = 1:length(zhist)
    z = zhist(i);
    C = Chist(i);
    S = Shist(i);
    if tol <= z || z <= - tol % Elliptical or Hyperbolic orbit
        dSdz = 0.5*(C - 3*S)/z; % Bate 5.3-20
    else % Parabolic orbit
        dSdz = - 1/factorial(5);
        for k = 1:K
            dSdz = dSdz ...
                + (k+1)*(-1)^(k-1)*z^k/factorial(2*k+5); % Bate 5.3-25
        end
    end
    dSdzhist(i) = dSdz;
end
end