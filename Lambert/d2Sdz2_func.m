function [ d2Sdz2hist ] = d2Sdz2_func( zhist, Chist, Shist, dCdzhist, ...
    dSdzhist, tol, K )
% Second Derivative of S parameter with respect to z
% Derivatives of the equation Bate 5.3-20 and 5.3-25 from "Fundamentals 
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
% dSdzhist - 1 x N double matrix
%            Second derivative of "Sine" parameter
% 
% @author: Matt Marti
% @date: 2019-03-28

% Preallocate
if numel(zhist) ~= 1
    d2Sdz2hist = zeros(size(zhist,1),size(zhist,2));
else
    d2Sdz2hist = 0;
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
        d2Sdz2 = 0.5*ooz*(dCdz - 3*dSdz) - 0.5*ooz*ooz*(C - 3*S);
    else % Parabolic orbit
        denom_i = 7*6*5*4*3*2;
        d2Sdz2 = 2/denom_i;
        for k = 1:K
            denom_i = denom_i * (2*k+7)*(2*k+6);
            d2Sdz2 = d2Sdz2 ...
                + (k+1)*(k+2)*(-z)^(k) / denom_i; % factorial(2*k+7);
        end
    end
    d2Sdz2hist(i) = d2Sdz2;
end
end