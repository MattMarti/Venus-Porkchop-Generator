function [v1vec, v2vec] = lambert(mu, r1vec, r2vec, tof, ccwflag, tol, K)% Lambert Algorithm to find velocity given time of flight% Computes velocity at both points of an orbit given the time of flight and% position vectors of the orbit.% % Based on the algorithm presented by Bate in "Fundamentals of% Astrodynamics", 1971.% % Also uses the DerAstrodynamics algorithm for determining whether the % transfer orbit has a change in true anomoly of greater than 180 degrees.% This function assumes planetary orbits are counter-clockwise% % @arg% mu        - double%             Gravitational Constant% r1vec     - 3 x 1 double matrix%             Start point Position vector% r2vec     - 3 x 1 double matrix%             Destination point position vector% tof       - double%             Time of flight% ccwflag   - bool (optional)%             Counter-clockwise orbit flag. True will cause the function to%             compute a counter-clockwise transfer orbit.%             Defaults to true if not specified.% tol       - double%             Toleranc of iterations% K         - int%             Maximum number of iterations% % @return% v1vec     - 3 x 1 double matrix%             Start point velocity vector% v2vec     - 3 x 1 double matrix%             Destination point velocity vector% % @author: Matt Marti% @date: 2019-03-26% Define square root of 2global PORKCHOPCODE_LAMBERT_SQRT2if isempty(PORKCHOPCODE_LAMBERT_SQRT2)    PORKCHOPCODE_LAMBERT_SQRT2 = sqrt(2);endsqrt2 = PORKCHOPCODE_LAMBERT_SQRT2;% Run test script if no argumentsif nargin == 0    Test_lambertendif nargin < 5    ccwflag = 1; % Defaultend% Constantsconstantsepsilon = 1e-2;% Check inputsassert(length(r1vec) == 3, 'incorrect size of position vector');assert(length(r2vec) == 3, 'incorrect size of position vector');assert(tof >= 0, 'Negative time of flight');% Compute normsr1 = norm(r1vec);r2 = norm(r2vec);% Compute geometry parameteralpha = dot([0;0;1], cross(r1vec,r2vec)); % DerAstrodynamics clockwise orbitcosdelnu = dot(r1vec,r2vec) / (r1*r2);delnu = acos(cosdelnu);DM = 1;if ccwflag % If counter-clockwise    if alpha < 0        delnu = 2*pi - delnu;        DM = -1;    endelse % Clockwise transfer    if alpha > 0        delnu = 2*pi - delnu;        DM = -1;    endendA = DM*sqrt(r1*r2*(1+cosdelnu)); % Bate 5.3-27tau = A / (r1 + r2);%% Initial guess% Not the real way to initialize the guess, but the algorithm is relatively% robustif hyperflag    k0 = 1.5 * sqrt2;else    k0 = 0 * sqrt2;end%% Lambert solution with Halley's Method Iterationi = 0;ki = k0;while i < IMAX_PITERATION        % Time of flight error    [tki, tpki, tppki] = t_func(ki, S, tau, epsilon);    Lki = tki - tof; % Eqn 31        % Lagrange parameters    kitau = ki*tau;    f = 1 - (1-kitau)/r1; % Eqn 32    g = tau * (r1+r2)*sqrt(1-kitau); % Eqn 33    gdot = 1 - (1-kitau) / r2; % Eqn 34        % Velocities    v1vec = (r2vec - f * r1vec) / g; % Eqn 35    v2vec = (gdot*r2vec - r1vec) / g; % Eqn 36        % End condition    if abs(Lki) < 1e-6        break;    end        % Increment k    deltak = -Lki / (tpki - 0.5*Lki*tppki/tpki); % Eqn 40    ki = ki + deltak;endend