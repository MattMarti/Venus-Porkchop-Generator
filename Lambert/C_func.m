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
if numel(zhist) ~= 1
    Chist = zeros(size(zhist,1),size(zhist,2));
else
    Chist = 0;
end

% Loop
for i = 1:length(zhist)
    z = zhist(i);
    if z <= - tol % Elliptical orbit
        Chist(i) = (1 - cos(sqrt(z))) ./ z;
    elseif tol <= z % Hyperbolic orbit
        Chist(i) = (1 - cosh(sqrt(-z))) ./ z;
    else % Series approximation for Parabolic orbit
        Chist(i) = 0;
        denom_i = 1;
        for k = 0:K
            denom_i = denom_i * (2*k+2) * (2*k+1);
            Chist(i) = Chist(i) + (-z)^k / denom_i; % factorial(2*k+2);
        end
    end
end
end