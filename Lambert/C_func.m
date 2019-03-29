function [ Chist ] = C_func( zhist, tol, K )
% "Cosine" parameter for Sundman Transform and the Lambert problem.
% Implementation of the equation Bate 4.4-10 from "Fundamentals of
% Astrodynamics", 1971.
% 
% @arg
% zhist - 1 x N double matrix
%         z parameter value
% tol   - double
%         Tolerance to apply the series approximation for parabolic orbit.
% K     - int
%         Series approximation number of iterations
% 
% @return
% Chist - 1 x N double matrix
%         "Cosine" parameter
% 
% @author: Matt Marti
% @date: 2019-03-28

% Preallocate
Chist = zeros(size(zhist));

% Loop
for i = 1:length(zhist)
    z = zhist(i);
    if z <= - tol % Elliptical orbit
        Chist(i) = (1 - cos(sqrt(z))) ./ z;
    elseif tol <= z % Hyperbolic orbit
        Chist(i) = (1 - cosh(sqrt(-z))) ./ z;
    else % Series approximation for Parabolic orbit
        Chist(i) = 0;
        for k = 0:K
            Chist(i) = Chist(i) + (-z).^k ./ factorial(2*k+2);
        end
    end
end
end