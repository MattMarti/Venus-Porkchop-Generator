function [xrv,dxrv_dxequinoctial] = equinoctial2rv(xequinoctial,mu)
%
%  Copyright (c) 2013, 2018 Mark L. Psiaki.  All rights reserved.  
%
%
%  This function transforms from the equinoctial coordinates defined
%  in the paper 
%
%    Broucke, R.A., and Cefloa, P.J., "On the Equinoctial Orbit Elements",
%    Celestial Mechanics, Vol. 5, 1972, pp. 303-310.
%  
%  to Cartesian position/velocity coordinates.  This function also
%  computes the Jacobian partial derivative of this transformation.
%  
%
%  Inputs:
%
%    xequinoctial        = [a;h;k;lambda;p;q], the 6-by-1 equinoctial
%                        satellite orbit state vector.  a is the semi-
%                        major axis, h = ecc*sin(omega+Omega),
%                        with ecc being the eccentricity, omega being
%                        the argument or perigee, and Omega being the
%                        longitude of the ascending node, k = ecc*...
%                        cos(omega+Omega), lambda = M + omega + Omega,
%                        with M being the mean anomly, p = ...
%                        tan(0.5*inc)*sin(Omega), with inc being the
%                        inclination, and q = tan(0.5*inc)*cos(Omega).
%                        Thus, xequinoctial(1,1) = a is in meters,
%                        xequinoctial([2,3,5,6],1) = [h;k;p;q]
%                        is non-dimensional, and xequinoctial(4,1) ...
%                        = lambda is given in units of radians.
%
%    mu                  The central gravitational constant of the body
%                        about which the satellite is orbiting, in
%                        units of m^3/sec^2.  mu = 3.986005e+14 m^3/sec^2
%                        is the value to use for orbits about the Earth.
%
%  Outputs:
%
%    xrv                 = [r;v], the 6-by-1 Cartesian position/velocity
%                        satellite orbit state vector.  r is the 3-by-1
%                        position vector in meters units, and v is the
%                        3-by-1 velocity vector in m/sec units.
%                        These two vectors are given in inertial
%                        ECIF coordinates.
%
%    dxrv_dxequinoctial  The 6-by-6 Jacobian partial derivative of xrv
%                        with respect to xequinoctial.  Its units are 
%                        consistent with those of xrv and xequinoctial.
%                        Thus, dxrv_dxequinoctial(1:3,1) is 
%                        dimensionless, dxrv_dxequinoctial(1:3,[2,3,5,6])
%                        has units of meters, dxrv_dxequinoctial(1:3,4)
%                        has units of meters/radian, 
%                        dxrv_dxequinoctial(4:6,1) has units of 1/sec,
%                        dxrv_dxequinoctial(4:6,[2,3,5,6]) has units of
%                        meters/sec, and dxrv_dxequinoctial(4:6,4)
%                        has units of meters/(radian-sec).
%

%
%  Set up and solve the equinoctial version of Kepler's equation for
%  F, the eccentric longitude, as a function of the mean longitude lambda
%  and of h and k.  Kepler's equation for F takes the form:
%
%       lambda = F + h*cos(F) - k*sin(F)
%
   h = xequinoctial(2,1);
   k = xequinoctial(3,1);
   lambda = xequinoctial(4,1);
   icount = 0;
   testdone0 = 0;
   F = lambda; % STUDENT RESPONSE
   facabsdelFlowlim = 1000*eps;
   absdelFlowlim = facabsdelFlowlim*2*pi;
   while testdone0 == 0
      icount = icount + 1;
      if icount > 150
         error(['Error in equinoctial2rv.m: Kepler solution',...
                ' failed to converge in 150 iterations.'])
      end
      sinF = sin(F);
      cosF = cos(F);
      g = lambda - F - h*cosF + k*sinF; % STUDENT RESPONSE
      dg_dF = -1 + h*cosF + k*cosF; % STUDENT RESPONSE
      deltaF = - g / dg_dF; % STUDENT RESPONSE
      F = F + deltaF;
      if abs(deltaF) <= max([absdelFlowlim;(facabsdelFlowlim*abs(F))])
         testdone0 = 1;
      end
   end
%
%  Compute X1 and Y1 after computing eccsq and various related factors.
%
   sinF = sin(F);
   cosF = cos(F);
   eccsq = h^2 + k^2; % STUDENT RESPONSE
   onemeccsq = 1 - eccsq;
   sqrtonemeccsq = sqrt(onemeccsq);
   oneplussqrtonemeccsq = 1 + sqrtonemeccsq;
   oneoveroneplussqrtonemeccsq = 1/oneplussqrtonemeccsq;
   lambdaminusF = h*cosF - k*sinF;
   faclambdaminusF = lambdaminusF*oneoveroneplussqrtonemeccsq;
   X1overa = cosF - k - h*faclambdaminusF; % STUDENT RESPONSE
   Y1overa = sinF - h + k*faclambdaminusF; % STUDENT RESPONSE
   a = xequinoctial(1,1);
   X1 = a*X1overa;
   Y1 = a*Y1overa;
%
%  Compute the transformation from orbital coordinates to the
%  final Cartesian coordinates.  Compute only the first 2
%  columns of it, the X1 column and the Y1 column, because
%  Z1 = 0 by definition.
%
   p = xequinoctial(5,1);
   q = xequinoctial(6,1);
   psq = p^2;
   qsq = q^2;
   psqplusqsq = psq + qsq;
   psqminusqsq = psq - qsq; 
   onepluspsqplusqsq = 1 + psqplusqsq;
   oneoveronepluspsqplusqsq = 1/onepluspsqplusqsq;
   twopq = 2*p*q;
   oneminuspsqplusqsq = 1 - psqminusqsq;
   onepluspsqminusqsq = 1 + psqminusqsq;
   X1hatECIFun = [oneminuspsqplusqsq;twopq;(-2*p)];
   X1hatECIF = oneoveronepluspsqplusqsq*X1hatECIFun;
   Y1hatECIFun = [twopq;onepluspsqminusqsq;(2*q)];
   Y1hatECIF = oneoveronepluspsqplusqsq*Y1hatECIFun;
%
%  Compute the final position vector.
%
   r = X1hatECIF * X1 + Y1hatECIF * Y1; % STUDENT RESPONSE
%
%  Determine lambdadot and Fdot.
%
   asq = a^2;
   acu = a*asq;
   oneoveracu = 1/acu;
   nsq = mu*oneoveracu;
   n = sqrt(nsq);
   lambdadot = n; % STUDENT RESPONSE
   oneminushsinFminuskcosF = 1 - h*sinF - k*cosF;
   oveoveroneminushsinFminuskcosF = 1/oneminushsinFminuskcosF;
   Fdot = lambdadot*oveoveroneminushsinFminuskcosF;
%
%  Compute X1dot and Y1dot based on lambdadot and Fdot.
%
   sinFdot = cosF*Fdot;
   cosFdot = - sinF*Fdot;
   lambdadotminusFdot = h*cosFdot - k*sinFdot;
   lambdadotminusFdotoveroneplussqrtonemeccsq = ...
             lambdadotminusFdot*oneoveroneplussqrtonemeccsq;
   X1overadot = -sinF*Fdot - h*lambdadotminusFdotoveroneplussqrtonemeccsq; % STUDENT RESPONSE
   Y1overadot = cosF*Fdot + k*lambdadotminusFdotoveroneplussqrtonemeccsq; % STUDENT RESPONSE
   X1dot = a*X1overadot;
   Y1dot = a*Y1overadot;
%
%  Compute the final velocity using X1dot, Y1dot, and the
%  transformation from orbital coordinates to the
%  final Cartesian coordinates.
%
   v = X1hatECIF * X1dot + Y1hatECIF * Y1dot; % STUDENT RESPONSE
%
%  Pack r and v into the xrv output.
%
   xrv = [r;v];
%
%  Compute the Jacobian partial derivative of r with 
%  respect to the equinoctial state.
%
%  First compute the partial derivative of F with respect ot the 
%  relevant elements of the equinoctial state.  It depends only on 
%  h, k, and lambda.
%
   dF_dh = oveoveroneminushsinFminuskcosF*(-cosF);
   dF_dk = oveoveroneminushsinFminuskcosF*(sinF);
   dF_dlambda = oveoveroneminushsinFminuskcosF;
%
%  Next compute the derivatives of X1 and Y1 with respect to the
%  relevant elements of the equinoctial state.  They depend only on 
%  a, h, k, and lambda.
%
   dsinF_dh = cosF*dF_dh;
   dsinF_dk = cosF*dF_dk;
   dsinF_dlambda = cosF*dF_dlambda;
   dcosF_dh = - sinF*dF_dh;
   dcosF_dk = - sinF*dF_dk;
   dcosF_dlambda = - sinF*dF_dlambda;
   deccsq_dh = 2*h;
   deccsq_dk = 2*k;
   donemeccsq_dh = - deccsq_dh;
   donemeccsq_dk = - deccsq_dk;
   halfosqrtonemeccsq = 0.5/sqrtonemeccsq;
   dsqrtonemeccsq_dh = halfosqrtonemeccsq*donemeccsq_dh;
   dsqrtonemeccsq_dk = halfosqrtonemeccsq*donemeccsq_dk;
   doneplussqrtonemeccsq_dh = dsqrtonemeccsq_dh;
   doneplussqrtonemeccsq_dk = dsqrtonemeccsq_dk;
   negoneooneoveroneplussqrtonemeccsqsq = ...
         - (oneoveroneplussqrtonemeccsq^2);
   doneoveroneplussqrtonemeccsq_dh = ...
         negoneooneoveroneplussqrtonemeccsqsq*doneplussqrtonemeccsq_dh;
   doneoveroneplussqrtonemeccsq_dk = ...
         negoneooneoveroneplussqrtonemeccsqsq*doneplussqrtonemeccsq_dk;
   dlambdaminusF_dh = cosF + h*dcosF_dh - k*dsinF_dh;
   dlambdaminusF_dk = h*dcosF_dk - sinF - k*dsinF_dk;
   dlambdaminusF_dlambda = h*dcosF_dlambda - k*dsinF_dlambda;
   dfaclambdaminusF_dh = ...
         dlambdaminusF_dh*oneoveroneplussqrtonemeccsq + ...
         lambdaminusF*doneoveroneplussqrtonemeccsq_dh;
   dfaclambdaminusF_dk = ...
         dlambdaminusF_dk*oneoveroneplussqrtonemeccsq + ...
         lambdaminusF*doneoveroneplussqrtonemeccsq_dk;
   dfaclambdaminusF_dlambda = ...
         dlambdaminusF_dlambda*oneoveroneplussqrtonemeccsq;
   dX1overa_dh = dcosF_dh - faclambdaminusF - h*dfaclambdaminusF_dh;
   dX1overa_dk = dcosF_dk - 1 - h*dfaclambdaminusF_dk;
   dX1overa_dlambda = dcosF_dlambda - h*dfaclambdaminusF_dlambda;
   dY1overa_dh = dsinF_dh - 1 + k*dfaclambdaminusF_dh;
   dY1overa_dk = dsinF_dk + faclambdaminusF + k*dfaclambdaminusF_dk;
   dY1overa_dlambda = dsinF_dlambda + k*dfaclambdaminusF_dlambda;
   dX1_da = X1overa;
   dX1_dh = a*dX1overa_dh;
   dX1_dk = a*dX1overa_dk;
   dX1_dlambda = a*dX1overa_dlambda;
   dY1_da = Y1overa;
   dY1_dh = a*dY1overa_dh;
   dY1_dk = a*dY1overa_dk;
   dY1_dlambda = a*dY1overa_dlambda;
%
%  Compute the partial derivatives with respect to p and q
%  of the transformation from orbital coordinates to the
%  final Cartesian coordinates.  Compute the partial 
%  derivatives only of the first 2 columns of it, the X1 
%  column and the Y1 column, because Z1 = 0 by definition.
%
   dpsq_dp = 2*p;
   dqsq_dq = 2*q;
   dpsqplusqsq_dp = dpsq_dp;
   dpsqplusqsq_dq = dqsq_dq;
   dpsqminusqsq_dp = dpsq_dp; 
   dpsqminusqsq_dq = - dqsq_dq;
   donepluspsqplusqsq_dp = dpsqplusqsq_dp;
   donepluspsqplusqsq_dq = dpsqplusqsq_dq;
   negoneoveronepluspsqplusqsqsq = - (oneoveronepluspsqplusqsq^2);
   doneoveronepluspsqplusqsq_dp = ...
         negoneoveronepluspsqplusqsqsq*donepluspsqplusqsq_dp;
   doneoveronepluspsqplusqsq_dq = ...
         negoneoveronepluspsqplusqsqsq*donepluspsqplusqsq_dq;
   dtwopq_dp = 2*q;
   dtwopq_dq = 2*p;
   doneminuspsqplusqsq_dp = - dpsqminusqsq_dp;
   doneminuspsqplusqsq_dq = - dpsqminusqsq_dq;
   donepluspsqminusqsq_dp = dpsqminusqsq_dp;
   donepluspsqminusqsq_dq = dpsqminusqsq_dq;
   dX1hatECIFun_dp = [doneminuspsqplusqsq_dp;dtwopq_dp;(-2)];
   dX1hatECIFun_dq = [doneminuspsqplusqsq_dq;dtwopq_dq;0];
   dX1hatECIF_dp = doneoveronepluspsqplusqsq_dp*X1hatECIFun + ...
                   oneoveronepluspsqplusqsq*dX1hatECIFun_dp;
   dX1hatECIF_dq = doneoveronepluspsqplusqsq_dq*X1hatECIFun + ...
                   oneoveronepluspsqplusqsq*dX1hatECIFun_dq;
   dY1hatECIFun_dp = [dtwopq_dp;donepluspsqminusqsq_dp;0];
   dY1hatECIFun_dq = [dtwopq_dq;donepluspsqminusqsq_dq;(2)];
   dY1hatECIF_dp = doneoveronepluspsqplusqsq_dp*Y1hatECIFun + ...
                   oneoveronepluspsqplusqsq*dY1hatECIFun_dp;  
   dY1hatECIF_dq = doneoveronepluspsqplusqsq_dq*Y1hatECIFun + ...
                   oneoveronepluspsqplusqsq*dY1hatECIFun_dq;  
%
%  Compute the Jacobian partial derivatives of the final position 
%  with respect to the equinoctial state and place them in the
%  appropriate elements of dxrv_dxequinoctial.
%
   dxrv_dxequinoctial = zeros(6,6);
   dxrv_dxequinoctial(1:3,:) = ...
                 [(dX1_da*X1hatECIF + dY1_da*Y1hatECIF),...
                  (dX1_dh*X1hatECIF + dY1_dh*Y1hatECIF),...
                  (dX1_dk*X1hatECIF + dY1_dk*Y1hatECIF),...
                  (dX1_dlambda*X1hatECIF + dY1_dlambda*Y1hatECIF),...
                  (X1*dX1hatECIF_dp + Y1*dY1hatECIF_dp),...
                  (X1*dX1hatECIF_dq + Y1*dY1hatECIF_dq)];
%
%  Compute the partial derivatives of lambdadot and Fdot
%  with respect to the relevant elements of the equinoctial
%  state.  They depend only on a, h, k, and lambda.
%
   doneoveracu_da = - 3*oneoveracu*(1/a);
   dnsq_da = mu*doneoveracu_da;
   halfon = 0.5/n;
   dn_da = halfon*dnsq_da;
   dlambdadot_da = dn_da;
   doneminushsinFminuskcosF_dh = - sinF - h*dsinF_dh - k*dcosF_dh;
   doneminushsinFminuskcosF_dk = - h*dsinF_dk - cosF - k*dcosF_dk;
   doneminushsinFminuskcosF_dlambda = - h*dsinF_dlambda - k*dcosF_dlambda;
   negoveoveroneminushsinFminuskcosFsq = ...
         - (oveoveroneminushsinFminuskcosF^2);
   doveoveroneminushsinFminuskcosF_dh = ...
         negoveoveroneminushsinFminuskcosFsq*doneminushsinFminuskcosF_dh;
   doveoveroneminushsinFminuskcosF_dk = ...
         negoveoveroneminushsinFminuskcosFsq*doneminushsinFminuskcosF_dk;
   doveoveroneminushsinFminuskcosF_dlambda = ...
         negoveoveroneminushsinFminuskcosFsq*...
             doneminushsinFminuskcosF_dlambda;
   dFdot_da = dlambdadot_da*oveoveroneminushsinFminuskcosF;
   dFdot_dh = lambdadot*doveoveroneminushsinFminuskcosF_dh;
   dFdot_dk = lambdadot*doveoveroneminushsinFminuskcosF_dk;
   dFdot_dlambda = lambdadot*doveoveroneminushsinFminuskcosF_dlambda;
%
%  Compute partial derivatives of X1dot and Y1dot with respect to 
%  the relevant elements of the equinoctial state.  They depend 
%  only on a, h, k, and lambda.
%
   dsinFdot_da = cosF*dFdot_da;
   dsinFdot_dh = dcosF_dh*Fdot + cosF*dFdot_dh;
   dsinFdot_dk = dcosF_dk*Fdot + cosF*dFdot_dk;
   dsinFdot_dlambda = dcosF_dlambda*Fdot + cosF*dFdot_dlambda;
   dcosFdot_da = - sinF*dFdot_da;
   dcosFdot_dh = - dsinF_dh*Fdot - sinF*dFdot_dh;
   dcosFdot_dk = - dsinF_dk*Fdot - sinF*dFdot_dk;
   dcosFdot_dlambda = - dsinF_dlambda*Fdot - sinF*dFdot_dlambda;
   dlambdadotminusFdot_da = h*dcosFdot_da - k*dsinFdot_da;
   dlambdadotminusFdot_dh = cosFdot + h*dcosFdot_dh - k*dsinFdot_dh;
   dlambdadotminusFdot_dk = h*dcosFdot_dk - sinFdot - k*dsinFdot_dk;
   dlambdadotminusFdot_dlambda = h*dcosFdot_dlambda - k*dsinFdot_dlambda;
   dlambdadotminusFdotoveroneplussqrtonemeccsq_da = ...
           dlambdadotminusFdot_da*oneoveroneplussqrtonemeccsq;
   dlambdadotminusFdotoveroneplussqrtonemeccsq_dh = ...
           dlambdadotminusFdot_dh*oneoveroneplussqrtonemeccsq + ...
           lambdadotminusFdot*doneoveroneplussqrtonemeccsq_dh;
   dlambdadotminusFdotoveroneplussqrtonemeccsq_dk = ...
           dlambdadotminusFdot_dk*oneoveroneplussqrtonemeccsq + ...
           lambdadotminusFdot*doneoveroneplussqrtonemeccsq_dk;
   dlambdadotminusFdotoveroneplussqrtonemeccsq_dlambda = ...
           dlambdadotminusFdot_dlambda*oneoveroneplussqrtonemeccsq;
   dX1overadot_da = dcosFdot_da - ...
                    h*dlambdadotminusFdotoveroneplussqrtonemeccsq_da;
   dX1overadot_dh = dcosFdot_dh - ...
                    lambdadotminusFdotoveroneplussqrtonemeccsq - ...
                    h*dlambdadotminusFdotoveroneplussqrtonemeccsq_dh;
   dX1overadot_dk = dcosFdot_dk - ...
                    h*dlambdadotminusFdotoveroneplussqrtonemeccsq_dk;
   dX1overadot_dlambda = dcosFdot_dlambda - ...
                    h*dlambdadotminusFdotoveroneplussqrtonemeccsq_dlambda;
   dY1overadot_da = dsinFdot_da + ...
                    k*dlambdadotminusFdotoveroneplussqrtonemeccsq_da;
   dY1overadot_dh = dsinFdot_dh + ...
                    k*dlambdadotminusFdotoveroneplussqrtonemeccsq_dh;
   dY1overadot_dk = dsinFdot_dk + ...
                    lambdadotminusFdotoveroneplussqrtonemeccsq + ...
                    k*dlambdadotminusFdotoveroneplussqrtonemeccsq_dk;
   dY1overadot_dlambda = dsinFdot_dlambda + ...
                    k*dlambdadotminusFdotoveroneplussqrtonemeccsq_dlambda;
   dX1dot_da = X1overadot + a*dX1overadot_da;
   dX1dot_dh = a*dX1overadot_dh;
   dX1dot_dk = a*dX1overadot_dk;
   dX1dot_dlambda = a*dX1overadot_dlambda;
   dY1dot_da = Y1overadot + a*dY1overadot_da;
   dY1dot_dh = a*dY1overadot_dh;
   dY1dot_dk = a*dY1overadot_dk;
   dY1dot_dlambda = a*dY1overadot_dlambda;
%
%  Compute the Jacobian partial derivatives of the final velocity 
%  with respect to the equinoctial state and place them in the
%  appropriate elements of dxrv_dxequinoctial.
%
   dxrv_dxequinoctial(4:6,:) = ...
                 [(dX1dot_da*X1hatECIF + dY1dot_da*Y1hatECIF),...
                  (dX1dot_dh*X1hatECIF + dY1dot_dh*Y1hatECIF),...
                  (dX1dot_dk*X1hatECIF + dY1dot_dk*Y1hatECIF),...
                  (dX1dot_dlambda*X1hatECIF + dY1dot_dlambda*Y1hatECIF),...
                  (X1dot*dX1hatECIF_dp + Y1dot*dY1hatECIF_dp),...
                  (X1dot*dX1hatECIF_dq + Y1dot*dY1hatECIF_dq)];