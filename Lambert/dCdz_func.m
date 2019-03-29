function [ dCdzhist ] = dCdz_func( zhist, Chist, Shist, tol, K )
% Derivative of C parameter with respect to z
% Implementation of the equation Bate 5.3-21 and 5.3-24 from "Fundamentals 
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
% dCdzhist - 1 x N double matrix
%            Derivative of "Cosine" parameter
% 
% @author: Matt Marti
% @date: 2019-03-28

% Preallocate
dCdzhist = zeros(size(zhist));

% Loop
for i = 1:length(zhist)
    z = zhist(i);
    C = Chist(i);
    S = Shist(i);
    if tol <= z || z <= - tol % Elliptical or Hyperbolic orbit
        dCdz = 0.5*(1 - z*S - 2*C)/z; % Bate 5.3-21
    else % Parabolic orbit
        dCdz = - 1/factorial(4);
        for k = 1:K
            dCdz = dCdz ...
                + (k+1)*(-1)^(k-1)*z^k/factorial(2*k+4); % Bate 5.3-24
        end
    end
    dCdzhist(i) = dCdz;
end
end