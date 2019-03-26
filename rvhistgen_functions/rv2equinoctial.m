function [xequinoctial,dxequinoctial_dxrv,n,dn_dxrv] = ...
                      rv2equinoctial(xrv,mu)
%
%  Copyright (c) 2013, 2017, 2018 Mark L. Psiaki.  All rights reserved.  
%
%
%  This function transforms to the equinoctial coordinates defined
%  in the paper 
%
%    Broucke, R.A., and Cefloa, P.J., "On the Equinoctial Orbit Elements",
%    Celestial Mechanics, Vol. 5, 1972, pp. 303-310.
%  
%  from Cartesian position/velocity coordinates.  This function also
%  computes the Jacobian partial derivative of this transformation.
%  
%
%  Inputs:
%
%    xrv                 = [r;v], the 6-by-1 Cartesian position/velocity
%                        satellite orbit state vector.  r is the 3-by-1
%                        position vector in meters units, and v is the
%                        3-by-1 velocity vector in m/sec units.
%                        These are both given in inertial ECIF
%                        coordinates.
%
%    mu                  The central gravitational constant of the body
%                        about which the satellite is orbiting, in
%                        units of m^3/sec^2.  mu = 3.986005e+14 m^3/sec^2
%                        is the value to use for orbits about the Earth.
%
%  Outputs:
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
%    dxequinoctial_dxrv  The 6-by-6 Jacobian partial derivative of 
%                        xequinoctial with respect to xrv.  Its units are 
%                        consistent with those of xrv and xequinoctial.
%                        Thus, dxequinoctial_dxrv(1,1:3) is 
%                        dimensionless, dxequinoctial_dxrv(1,4:6)
%                        has units of seconds, dxequinoctial_dxrv(...
%                        [2,3,5,6],1:3) has units of 1/meter, 
%                        dxequinoctial_dxrv([2,3,5,6],4:6) has units of 
%                        sec/meter, dxequinoctial_dxrv(4,1:3) has
%                        units of radians/meter, and 
%                        dxequinoctial_dxrv(4,4:6) has units of 
%                        (radians-sec)/meter.
%
%    n                   The mean motion in rad/sec.
%
%    dn_dxrv             The 1-by-6 partial derivative of n with
%                        respect to xrv.  dn_dxrv(1,1:3) is
%                        given in units of rad/(sec-m), and
%                        dn_dxrv(1,4:6) is given in units of rad/m.
%

%
%  Determine the geocentric radius and the square of the
%  the velocity and use the vis viva equation in order to compute a.
%
   r = xrv(1:3,1);
   rmagsq = sum(r.^2);
   rmag = sqrt(rmagsq);
   v = xrv(4:6,1);
   vmagsq = sum(v.^2);
   oneoverrmag = 1/rmag;
   twooverrmag = 2*oneoverrmag;
   oneovermu = 1/mu;
   vmagsqovermu = vmagsq*oneovermu;
   oneovera = twooverrmag - vmagsqovermu; % STUDENT RESPONSE
   a = 1/oneovera;
%
%  Determine the angular-momentum-per-unit-mass vector
%  and its corresponding unit direction vector.
%
   hAM = cross(r,v); % STUDENT RESPONSE
   hAMmagsq = sum(hAM.^2);
   hAMmag = sqrt(hAMmagsq);
   oneoverhAMmag = 1/hAMmag;
   hAMhat = hAM*oneoverhAMmag;
%
%  Use the angular momemtum unit direction vector in order to 
%  determine p and q.
%
   oneminushAMhat3 = 1 - hAMhat(3,1);
   oneplushAMhat3 = 1 + hAMhat(3,1);
   oneoveroneplushAMhat3 = 1/oneplushAMhat3;
   psqplusqsq = oneminushAMhat3*oneoveroneplushAMhat3;
   oneplusppsqplusqsq = 1 + psqplusqsq;
   half_oneplusppsqplusqsq = 0.5*oneplusppsqplusqsq;
   p = hAMhat(1,1)*half_oneplusppsqplusqsq; % STUDENT RESPONSE
   q = -hAMhat(2,1)*half_oneplusppsqplusqsq; % STUDENT RESPONSE
%
%  Compute the transformation from orbital coordinates to the
%  final Cartesian coordinates.  Compute only the first 2
%  columns of it, the X1 column and the Y1 column, because
%  Z1 = 0 by definition.
%
   psq = p^2;
   qsq = q^2;
   psqminusqsq = psq - qsq; 
   oneoveroneplusppsqplusqsq = 1/oneplusppsqplusqsq;
   twopq = 2*p*q;
   oneminuspsqplusqsq = 1 - psqminusqsq;
   onepluspsqminusqsq = 1 + psqminusqsq;
   X1hatECIFun = [oneminuspsqplusqsq;twopq;(-2*p)];
   X1hatECIF = oneoveroneplusppsqplusqsq*X1hatECIFun;
   Y1hatECIFun = [twopq;onepluspsqminusqsq;(2*q)];
   Y1hatECIF = oneoveroneplusppsqplusqsq*Y1hatECIFun;   
%
%  Compute X1dot and Y1dot.
%
   X1dot = X1hatECIF' * v; % STUDENT RESPONSE
   Y1dot = Y1hatECIF' * v; % STUDENT RESPONSE
%
%  Compute n, the orbital rate.
%
   asq = a^2;
   acu = a*asq;
   oneoveracu = 1/acu;
   nsq = mu*oneoveracu;
   n = sqrt(nsq);
%
%  Compute the quantities sqrt(1 - ecc^2) and (1 + sqrt(1 - ecc^2)).
%
   Energy = -mu/(2*a); % STUDENT RESPONSE
   twoovermusq = 2*(oneovermu^2);
   Energy_hAMmagsq = Energy*hAMmagsq;
   eccsqminus1 = twoovermusq*Energy_hAMmagsq;
   oneminuseccsq = - eccsqminus1;
   sqrtoneminuseccsq = sqrt(oneminuseccsq);
   oneplussqrtoneminuseccsq = 1 + sqrtoneminuseccsq;
%
%  Compute the time rate of change of rmag.
%
   rdotv = sum(r.*v);
   rmagdot = rdotv*oneoverrmag;
%
%  Compute h and k by setting up and solving
%  a system of 2 linear equation in two unknowns.
%
%  It takes the form:
%
%      cospsi*h + sinpsi*k = u
%
%     -sinpsi*h + cospsi*k = w
%
%  where
%
%      cospsi = Y1dot/vmag
%
%      sinpsi = X1dot/vmag
%
%      u = - rmagdot/vmag
%
%      w = ((a/rmag) - 1)*sqrt(1 - ecc^2)*((n*a)/vmag)
%
%  The 2-by-2 matrix [cospsi,sinpsi;-sinpsi,cospsi]
%  is orthonormal.  Therefore, its inverse equals its
%  transpose, and the solution of this system of equations
%  is easy to compute.  It is easy to prove that the requirement
%  cospsi^2 + sinpsi^2 = 1 is met because vmag^2 = X1dot^2 + Y1dot^2.
%
   vmag = sqrt(vmagsq);
   oneovervmag = 1/vmag;
   cospsi = Y1dot/vmag; % STUDENT RESPONSE
   sinpsi = X1dot/vmag; % STUDENT RESPONSE
   u = -rmagdot*oneovervmag;
   aoverrmag = a*oneoverrmag;
   aoverrmagminusone = aoverrmag - 1;
   facw = aoverrmagminusone*sqrtoneminuseccsq;
   na = n*a;
   naovervmag = na*oneovervmag;
   w = facw*naovervmag;
   h = cospsi*u - sinpsi*w; % STUDENT RESPONSE
   k = sinpsi*u + cospsi*w; % STUDENT RESPONSE
%
%  Finish computing F.  Note that sinFfac and cosFfac
%  are different from, respectively, sinF and cosF through
%  multiplication by a common positive scalar.  It is not 
%  necessary to know or determine that scalar in order
%  to use sinFfac and cosFfac to compute F.  It might be 
%  helpful to one's understanding of this function, however,
%  to determine theoretically what that scalar is.  
%
   fachksinFcosF = aoverrmagminusone*na;
   sinFfac = h*fachksinFcosF - oneplussqrtoneminuseccsq*X1dot;
   cosFfac = k*fachksinFcosF + oneplussqrtoneminuseccsq*Y1dot;
   F = atan2(sinFfac, cosFfac); % STUDENT RESPONSE
%
%  Compute lambda.  Note that the method used below to
%  compute sinF and cosF could have been replaced by
%  sinF = sin(F) and cosF = cos(F) because F has already 
%  been calculated at this point of the code.  The method 
%  used below has been chosen because it is thought to be 
%  more efficient computationally.
%
   sinFfacsqpluscosFfacsq = sinFfac^2 + cosFfac^2;
   sqrtsinFfacsqpluscosFfacsq = sqrt(sinFfacsqpluscosFfacsq);
   oneoversqrtsinFfacsqpluscosFfacsq = 1/sqrtsinFfacsqpluscosFfacsq;
   cosF = cosFfac*oneoversqrtsinFfacsqpluscosFfacsq;
   sinF = sinFfac*oneoversqrtsinFfacsqpluscosFfacsq;
   lambda = F + h*cosF - k*sinF; % STUDENT RESPONSE
%
%  Pack the equinoctial states into the output vector xequinoctial.
%
   xequinoctial = [a;h;k;lambda;p;q];
%
%  Now do the derivatives with respect to xrv.  Efficiency
%  is not the most important issue pursued in this derivation.
%  Instead, it is designed to parallel the developments given above.
%  Therefore, some derivatives with respect to xrv are carried out
%  with either the first 3 columns or the last 3 columns being
%  all zeros due to the quantity in question depending only
%  on, respectively, the velocity v or the positiion r.  As the
%  calculations progress, there are fewer situations in which
%  this is the case.
%
   dr_dxrv = [eye(3),zeros(3,3)];
   rtr = r';
   drmagsq_dxrv = (2*rtr)*dr_dxrv;
   halfoverrmag = 0.5*oneoverrmag;
   drmag_dxrv = halfoverrmag*drmagsq_dxrv;
   dv_dxrv = [zeros(3,3),eye(3)];
   vtr = v';
   dvmagsq_dxrv = (2*vtr)*dv_dxrv;
   negoneoverrmagsq = - (oneoverrmag^2);
   doneoverrmag_dxrv = negoneoverrmagsq*drmag_dxrv;
   dtwooverrmag_dxrv = 2*doneoverrmag_dxrv;
   dvmagsqovermu_dxrv = dvmagsq_dxrv*oneovermu;
   doneovera_dxrv = dtwooverrmag_dxrv - dvmagsqovermu_dxrv;
   negasq = - asq;
   da_dxrv = negasq*doneovera_dxrv;
%
   rCPEmat = [      0, -r(3,1),  r(2,1);...
               r(3,1),       0, -r(1,1);...
              -r(2,1),  r(1,1),       0];
   vCPEmat = [      0, -v(3,1),  v(2,1);...
               v(3,1),       0, -v(1,1);...
              -v(2,1),  v(1,1),       0];
   dhAM_dxrv = rCPEmat*dv_dxrv - vCPEmat*dr_dxrv;
   dhAMmagsq_dxrv = (2*(hAM'))*dhAM_dxrv;
   halfoverhAMmag = 0.5*oneoverhAMmag;
   dhAMmag_dxrv = halfoverhAMmag*dhAMmagsq_dxrv;
   negoneoverhAMmagsq = - (oneoverhAMmag^2);
   doneoverhAMmag_dxrv = negoneoverhAMmagsq*dhAMmag_dxrv;
   dhAMhat_dxrv = dhAM_dxrv*oneoverhAMmag + ...
                  hAM*doneoverhAMmag_dxrv;
%
   doneminushAMhat3_dxrv = - dhAMhat_dxrv(3,:);
   doneplushAMhat3_dxrv = dhAMhat_dxrv(3,:);
   negoneoveroneplushAMhat3sq = - (oneoveroneplushAMhat3^2);
   doneoveroneplushAMhat3_dxrv = ...
           negoneoveroneplushAMhat3sq*doneplushAMhat3_dxrv;
   dpsqplusqsq_dxrv = doneminushAMhat3_dxrv*oneoveroneplushAMhat3 + ...
                      oneminushAMhat3*doneoveroneplushAMhat3_dxrv;
   doneplusppsqplusqsq_dxrv = dpsqplusqsq_dxrv;
   dhalf_oneplusppsqplusqsq_dxrv = 0.5*doneplusppsqplusqsq_dxrv;
   dp_dxrv = dhAMhat_dxrv(1,:)*half_oneplusppsqplusqsq + ...
             hAMhat(1,1)*dhalf_oneplusppsqplusqsq_dxrv;
   dq_dxrv = - dhAMhat_dxrv(2,:)*half_oneplusppsqplusqsq - ...
               hAMhat(2,1)*dhalf_oneplusppsqplusqsq_dxrv;
%
   dpsq_dp = 2*p;
   dqsq_dq = 2*q;
   dpsqplusqsq_dp = dpsq_dp;
   dpsqplusqsq_dq = dqsq_dq;
   dpsqminusqsq_dp = dpsq_dp; 
   dpsqminusqsq_dq = - dqsq_dq;
   doneplusppsqplusqsq_dp = dpsqplusqsq_dp;
   doneplusppsqplusqsq_dq = dpsqplusqsq_dq;
   negoneoveroneplusppsqplusqsqsq = - (oneoveroneplusppsqplusqsq^2);
   doneoveroneplusppsqplusqsq_dp = ...
            negoneoveroneplusppsqplusqsqsq*doneplusppsqplusqsq_dp;
   doneoveroneplusppsqplusqsq_dq = ...
           negoneoveroneplusppsqplusqsqsq*doneplusppsqplusqsq_dq;
   dtwopq_dp = 2*q;
   dtwopq_dq = 2*p;
   doneminuspsqplusqsq_dp = - dpsqminusqsq_dp;
   doneminuspsqplusqsq_dq = - dpsqminusqsq_dq;
   donepluspsqminusqsq_dp = dpsqminusqsq_dp;
   donepluspsqminusqsq_dq = dpsqminusqsq_dq;
   dX1hatECIFun_dp = [doneminuspsqplusqsq_dp;dtwopq_dp;(-2)];
   dX1hatECIFun_dq = [doneminuspsqplusqsq_dq;dtwopq_dq;0];
   dX1hatECIF_dp = doneoveroneplusppsqplusqsq_dp*X1hatECIFun + ...
                   oneoveroneplusppsqplusqsq*dX1hatECIFun_dp;
   dX1hatECIF_dq = doneoveroneplusppsqplusqsq_dq*X1hatECIFun + ...
                   oneoveroneplusppsqplusqsq*dX1hatECIFun_dq;
   dY1hatECIFun_dp = [dtwopq_dp;donepluspsqminusqsq_dp;0];
   dY1hatECIFun_dq = [dtwopq_dq;donepluspsqminusqsq_dq;(2)];
   dY1hatECIF_dp = doneoveroneplusppsqplusqsq_dp*Y1hatECIFun + ...
                   oneoveroneplusppsqplusqsq*dY1hatECIFun_dp;  
   dY1hatECIF_dq = doneoveroneplusppsqplusqsq_dq*Y1hatECIFun + ...
                   oneoveroneplusppsqplusqsq*dY1hatECIFun_dq; 
   dX1hatECIF_dxrv = dX1hatECIF_dp*dp_dxrv + dX1hatECIF_dq*dq_dxrv;
   dY1hatECIF_dxrv = dY1hatECIF_dp*dp_dxrv + dY1hatECIF_dq*dq_dxrv;
%
   X1hatECIFtr = X1hatECIF';
   Y1hatECIFtr = Y1hatECIF';
   dX1dot_dxrv = X1hatECIFtr*dv_dxrv + vtr*dX1hatECIF_dxrv;
   dY1dot_dxrv = Y1hatECIFtr*dv_dxrv + vtr*dY1hatECIF_dxrv;
%
   doneoveracu_dxrv = (-3*oneoveracu*oneovera)*da_dxrv;
   dnsq_dxrv = mu*doneoveracu_dxrv;
   halfovern = 0.5/n;
   dn_dxrv = halfovern*dnsq_dxrv;
%
   dEnergy_dxrv = 0.5*dvmagsq_dxrv - mu*doneoverrmag_dxrv;
   dEnergy_hAMmagsq_dxrv = dEnergy_dxrv*hAMmagsq + ...
                           Energy*dhAMmagsq_dxrv;
   deccsqminus1_dxrv = twoovermusq*dEnergy_hAMmagsq_dxrv;  
   doneminuseccsq_dxrv = - deccsqminus1_dxrv;
   halfoversqrtoneminuseccsq = 0.5/sqrtoneminuseccsq;
   dsqrtoneminuseccsq_dxrv = ...
              halfoversqrtoneminuseccsq*doneminuseccsq_dxrv;
   doneplussqrtoneminuseccsq_dxrv = dsqrtoneminuseccsq_dxrv;
%
   drdotv_dxrv = rtr*dv_dxrv + vtr*dr_dxrv;
   drmagdot_dxrv = drdotv_dxrv*oneoverrmag + rdotv*doneoverrmag_dxrv;
%
   halfovervmag = 0.5*oneovervmag;
   dvmag_dxrv = halfovervmag*dvmagsq_dxrv;
   negoneovervmagsq = - (oneovervmag^2);
   doneovervmag_dxrv = negoneovervmagsq*dvmag_dxrv;
   dcospsi_dxrv = dY1dot_dxrv*oneovervmag + Y1dot*doneovervmag_dxrv;
   dsinpsi_dxrv = dX1dot_dxrv*oneovervmag + X1dot*doneovervmag_dxrv;
   du_dxrv = - drmagdot_dxrv*oneovervmag - rmagdot*doneovervmag_dxrv;
   daoverrmag_dxrv = da_dxrv*oneoverrmag + a*doneoverrmag_dxrv;
   daoverrmagminusone_dxrv = daoverrmag_dxrv;
   dfacw_dxrv = daoverrmagminusone_dxrv*sqrtoneminuseccsq + ...
                aoverrmagminusone*dsqrtoneminuseccsq_dxrv;
   dna_dxrv = dn_dxrv*a + n*da_dxrv;
   dnaovervmag_dxrv = dna_dxrv*oneovervmag + na*doneovervmag_dxrv;
   dw_dxrv = dfacw_dxrv*naovervmag + facw*dnaovervmag_dxrv;
   dh_dxrv = dcospsi_dxrv*u + cospsi*du_dxrv - dsinpsi_dxrv*w - ...
             sinpsi*dw_dxrv;
   dk_dxrv = dsinpsi_dxrv*u + sinpsi*du_dxrv + dcospsi_dxrv*w + ...
             cospsi*dw_dxrv;  
%
   dfachksinFcosF_dxrv = daoverrmagminusone_dxrv*na + ...
                         aoverrmagminusone*dna_dxrv;
   dsinFfac_dxrv = dh_dxrv*fachksinFcosF + h*dfachksinFcosF_dxrv - ...
                   doneplussqrtoneminuseccsq_dxrv*X1dot - ...
                   oneplussqrtoneminuseccsq*dX1dot_dxrv;                     
   dcosFfac_dxrv = dk_dxrv*fachksinFcosF + k*dfachksinFcosF_dxrv + ...
                   doneplussqrtoneminuseccsq_dxrv*Y1dot + ...
                   oneplussqrtoneminuseccsq*dY1dot_dxrv;
   oneoversinFfacsqpluscosFfacsq = 1/sinFfacsqpluscosFfacsq;
   dF_dxrv = (-sinFfac*dcosFfac_dxrv + cosFfac*dsinFfac_dxrv)*...
              oneoversinFfacsqpluscosFfacsq;
%
   dcosF_dxrv = - sinF*dF_dxrv;
   dsinF_dxrv = cosF*dF_dxrv;
   dlambda_dxrv = dF_dxrv + dh_dxrv*cosF + h*dcosF_dxrv - ...
                  dk_dxrv*sinF - k*dsinF_dxrv;
%
%  Pack the equinoctial state partial derivatives into the output 
%  array dxequinoctial_dxrv.
%
   dxequinoctial_dxrv = [da_dxrv;dh_dxrv;dk_dxrv;dlambda_dxrv;...
                         dp_dxrv;dq_dxrv];