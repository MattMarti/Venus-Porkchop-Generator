function [xrvhist,dxrvhist_dxrvepoch,dxrvhist_dthist] = ...
                      rvhistgen(thist,tepoch,xrvepoch,mu)
%
%  Copyright (c) 2017, 2018 Mark L. Psiaki.  All rights reserved.  
%
%
%  This function computes the inertial/ECIF positions and 
%  velocities along a Keplerian orbit at a set of input times.
%  It does so by using non-singular equinoctial orbital
%  elements as defined in the paper.
%
%    Broucke, R.A., and Cefloa, P.J., "On the Equinoctial Orbit Elements",
%    Celestial Mechanics, Vol. 5, 1972, pp. 303-310.
%  
%  It first converts from position and velocity at an epoch time
%  to equinoctial elements.  It then propagates the equinoctial
%  mean longitude to the output times.  Finally, it transforms
%  the resulting equinoctial elements at the output times back
%  into ECIF position and velocity vectors.  This function uses
%  the function rv2equinoctial.m in order to transform from
%  ECIF position and velocity into equinoctial elements.  It
%  uses the function equinoctial2rv.m in order to transform
%  from equinoctial elements back into ECIF position and velocity.
%
%  This functon also computes the differential sensitivities of
%  the orbital position/velocity state time histories to the 
%  changes in the input position/velocity state elements at the 
%  epoch time.  It also computes the differential sensitivities 
%  of the position/velocity state time histories to changes in
%  their output times.
%  
%  
%
%  Inputs:
%
%    thist               = [t1;t2;t3;...;tNsamps], the Nsamps-by-1
%                        vector of times at which the position/
%                        velocity state vector will be computed
%                        and output in xrvhist.  These times are
%                        given in seconds measured relative to the 
%                        start of a particular GPS week.
%
%    tepoch              The scalar epoch time at which the position/
%                        velocity state in xrvepoch is defined.  This 
%                        time is given in seconds measured relative to 
%                        the start of a particular GPS week.
%
%    xrvepoch            = [repoch;vepoch], the 6-by-1 Cartesian 
%                        position/velocity satellite orbit state 
%                        vector that applies at time tepoch.  repoch
%                        is the 3-by-1 position vector in meters units,
%                        and vepoch is the 3-by-1 velocity vector in 
%                        meters/second units.  These vectors are  
%                        both given in inertial ECIF coordinates.
%
%    mu                  The central gravitational constant of the body
%                        about which the satellite is orbiting, in
%                        units of m^3/sec^2.  mu = 3.986005e+14 m^3/sec^2
%                        is the value to use for orbits about the Earth.
%
%  Outputs:
%
%    xrvhist             = [xrv1';xrv2';xrv3';...;xrvNsamps'], the
%                        Nsamps-by-6 array of ECIF satellite 
%                        position/velocity vectors at the times
%                        in thist.  xrvk = xrvhist(k,:)' = ...
%                        [rk;vk] is the 6-by-1 position/velocity
%                        state that applies at time tk = thist(k,1),
%                        with rk given in meters units and vk 
%                        given in meters/second units.  Thus,
%                        xrvhist(:,1) is the Nsamps-by-1 time
%                        history of the ECIF X component of the
%                        satellite position, xrvhist(:,2) is the
%                        ECIF Y position component time history,
%                        xrvhist(:,3) is the ECIF Z position  
%                        component time history, xrvhist(:,4) is 
%                        the ECIF X velocity component time history,
%                        xrvhist(:,5) is the ECIF Y velocity  
%                        component time history, and xrvhist(:,6) is 
%                        the ECIF Z velocity component time history.
%
%    dxrvhist_dxrvepoch  The Nsamps-by-6-by-6 array that gives
%                        the first partial derivatives of 
%                        xrvhist with respect to the components
%                        of xrvepoch.  The Nsamps-by-6 array
%                        dxrvhist_dxrvepoch(:,:,jj) is the first
%                        partial derivative of xrvhist with
%                        respect to xrvepoch(jj,1).  This
%                        is true for all jj = 1:6.  Thus,
%                        dxrvhist_dxrvepoch(:,1:3,jj) is
%                        non-dimensional for all jj = 1:3, and
%                        dxrvhist_dxrvepoch(:,4:6,jj) has
%                        units of 1/seconds for all jj = 1:3.
%                        Similarly, dxrvhist_dxrvepoch(:,1:3,jj)
%                        has units of seconds for all jj = 4:6,
%                        and dxrvhist_dxrvepoch(:,4:6,jj) is
%                        non-dimensional for all jj = 4:6.
%
%    dxrvhist_dthist     The Nsamps-by-6 array that gives the
%                        first partial derivatives of the elements 
%                        of xrvhist with respect to changes in
%                        their corresponding output times
%                        in thist.  Thus, dxrvk_dtk = ...
%                        dxrvhist_dthist(k,:)' is the 6-by-1
%                        first partial derivative of xrvk with
%                        respect to tk.  The units of
%                        dxrvhist_dthist(:,1:3) are meters/second,
%                        and the units of dxrvhist_dthist(:,4:6)
%                        are meters/second^2.  One might normally
%                        consider these derivatives to be total
%                        derivatives rather than partial derivatives.
%                        The language "partial derivative" is
%                        used here in order to specify that the
%                        input xrvepoch is held constant while
%                        the given value in thist is considered to
%                        be variable in each derivative calculation.
%

%
%  Determine the equinoctial elements at the epoch time.
%
   [xequinoctialepoch,dxequinoctialepoch_dxrvepoch,n,dn_dxrvepoch] = ...
          rv2equinoctial(xrvepoch,mu); % STUDENT RESPONSE
%
%  Determine the number of output sample times.
%
   Nsamps = size(thist,1);
%
%  Set up output arrays that are initially populated by
%  dummy zero values.
%
   xrvhist = zeros(Nsamps,6);
   dxrvhist_dxrvepoch = zeros(Nsamps,6,6);
   dxrvhist_dthist = zeros(Nsamps,6);
%
%  Retrieve the mean longitude at epoch and its partial
%  derivative with respect to xrvepoch.
%
   lambdaepoch = xequinoctialepoch(4); % STUDENT RESPONSE
   dlambdaepoch_dxrvepoch = dxequinoctialepoch_dxrvepoch(4,:);
%
%  Work through the output sample times in this while loop, one
%  time per iteration of the loop.  Compute the new equinoctial
%  elements at each time and then transform them back to
%  position/velocity components.  Store the results in the
%  appropriate row of xrvhist.  Afterward, calculate the needed
%  partial derivatives and store them in the relevant rows of 
%  dxrvhist_dxrvepoch and dxrvhist_dthist.
%
   for k = 1:Nsamps
      tk = thist(k,1);
      deltk = tk - tepoch;
      lambdak = lambdaepoch + deltk*n; % STUDENT RESPONSE
      xequinoctialk = [xequinoctialepoch(1:3,1);lambdak;...
                       xequinoctialepoch(5:6,1)]; % STUDENT RESPONSE
      [xrvk,dxrvk_dxequinoctialk] = equinoctial2rv(xequinoctialk,mu); % STUDENT RESPONSE
      xrvhist(k,:) = xrvk';
%
%  Calculate partial derivatives.
%
      dlambdak_dxrvepoch = dlambdaepoch_dxrvepoch + dn_dxrvepoch*deltk;
      dxequinoctialk_dxrvepoch = ...
                    [dxequinoctialepoch_dxrvepoch(1:3,:);...
                     dlambdak_dxrvepoch;...
                     dxequinoctialepoch_dxrvepoch(5:6,:)];
      dxrvk_dxrvepoch = dxrvk_dxequinoctialk*dxequinoctialk_dxrvepoch;
      dxrvhist_dxrvepoch(k,:,:) = dxrvk_dxrvepoch;
%
      dlambdak_dtk = n;
      dxrvk_dlambdak = dxrvk_dxequinoctialk(:,4);
      dxrvk_dtk = dxrvk_dlambdak*dlambdak_dtk;
      dxrvhist_dthist(k,:) = dxrvk_dtk';
   end