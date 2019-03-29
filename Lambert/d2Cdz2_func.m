function [ d2Cdz2hist ] = d2Cdz2_func( zhist, Chist, Shist, dCdzhist, ...
    dSdzhist, tol, K )
% Derivative of C parameter with respect to z
% Derivatives of the equation Bate 5.3-21 and 5.3-24 from "Fundamentals 
% of Astrodynamics", 1971.
% 
% @arg
% zhist    - 1 x N double matrix
%            z parameter value
% Chist    - 1 x N double matrix
%            C parameter value
% Shist    - 1 x N double matrix
%            S parameter value
% dCdzhist - 1 x N double matrix
%            Derivative of C parameter value with respect to z
% dSdzhist - 1 x N double matrix
%            Derivative of S parameter value with respect to z
% tol      - double
%            Tolerance to apply the series approximation for parabolic orbit.
% K        - int
%            Series approximation number of iterations
% 
% @return
% dCdzhist - 1 x N double matrix
%            Second derivative of "Cosine" parameter
% 
% @author: Matt Marti
% @date: 2019-03-28

% Preallocate
if numel(zhist) ~= 1
    d2Cdz2hist = zeros(size(zhist,1),size(zhist,2));
else
    d2Cdz2hist = 0;
end

% Loop
for i = 1:length(zhist)
    z = zhist(i);
    C = Chist(i);
    S = Shist(i);
    dCdz = dCdzhist(i);
    dSdz = dSdzhist(i);
    if tol <= z || z <= - tol % Elliptical or Hyperbolic orbit
        ooz = 1/z;
        d2Cdz2 = 0.5*ooz*ooz*(z*S + 2*C - 1) ...
            - 0.5*ooz*(S + z*dSdz + 2*dCdz);
    else % Parabolic orbit
        denom_i = 6*5*4*3*2;
        d2Cdz2 = 2/denom_i;
        for k = 1:K
            denom_i = denom_i * (2*k+6)*(2*k+5);
            d2Cdz2 = d2Cdz2 ...
                + (k+1)*(k+2)*(-z)^k / denom_i; % factorial(2*k+6);
        end
    end
    d2Cdz2hist(i) = d2Cdz2;
end
end