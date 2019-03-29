function [ Shist ] = S_func( zhist, tol, K )
% "Sine" parameter for Sundman Transform and the Lambert problem.
% Implementation of the equation Bate 4.4-11 from "Fundamentals of
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
% Shist - 1 x N double matrix
%         "Sine" parameter
% 
% @author: Matt Marti
% @date: 2019-03-28

% Preallocate
Shist = zeros(size(zhist));

% Loop
for i = 1:length(zhist)
    z = zhist(i);
    if z <= - tol % Elliptical orbit
        sqrtz = sqrt(z);
        Shist(i) = (sqrtz - sin(sqrtz))./(sqrtz.^3);
    elseif tol <= z % Hyperbolic orbit
        sqrtmz = sqrt(-z);
        Shist(i) = (sinh(sqrtmz) - sqrtmz)./(sqrtmz.^3);
    else % Series approximation for Parabolic orbit
        Shist(i) = 0;
        for k = 0:K
            Shist(i) = Shist(i) + (-z).^k ./ factorial(2*k+3);
        end
    end
end


end