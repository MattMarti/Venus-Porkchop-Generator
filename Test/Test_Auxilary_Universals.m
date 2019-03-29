%% Test_Auxilary_Universals
% 
% Tests the "dCdz_func.m", "d2Cdz2_func.m", "dSdz_func.m", and 
% "d2Sdz2_func.m" to make sure that the derivatives are calculated 
% correctly
% 
% @author: Matt Marti
% @date: 2019-03-28

clear


%% Test 1: Elliptical Case

% Choose z > 0 for ellipse
z = 0.5*pi;
dz = 1e-4;
tol = 1e-2;
K = 16;

% Function at z
C = C_func(z, tol, K);
S = S_func(z, tol, K);

% Numerical 1st derivative
C_p = C_func(z+dz, tol, K);
C_m = C_func(z-dz, tol, K);
S_p = S_func(z+dz, tol, K);
S_m = S_func(z-dz, tol, K);
dCdz_true = 0.5*(C_p - C_m)/dz;
dSdz_true = 0.5*(S_p - S_m)/dz;

% Analytical Derivative
dCdz = dCdz_func(z, C, S, tol, K);
dSdz = dSdz_func(z, C, S, tol, K);

% Compare analytical and numerical derivative
assert(abs(dCdz_true - dCdz) < 1e-6, 'Bad derivative');
assert(abs(dSdz_true - dSdz) < 1e-6, 'Bad derivative');

% Numercial 2nd derivative
dCdz_p = dCdz_func(z+dz, C_p, S_p, tol, K);
dCdz_m = dCdz_func(z-dz, C_m, S_m, tol, K);
dSdz_p = dSdz_func(z+dz, C_p, S_p, tol, K);
dSdz_m = dSdz_func(z-dz, C_m, S_m, tol, K);
d2Cdz2_true = 0.5*(dCdz_p - dCdz_m)/dz;
d2Sdz2_true = 0.5*(dSdz_p - dSdz_m)/dz;

% Analytical 2nd derivative
d2Cdz2 = d2Cdz2_func(z, C, S, dCdz, dSdz, tol, K);
d2Sdz2 = d2Sdz2_func(z, C, S, dCdz, dSdz, tol, K);

% Compare analytical and numerical second derivative
assert(abs(d2Cdz2_true - d2Cdz2) < 1e-6, 'Bad derivative');
assert(abs(d2Sdz2_true - d2Sdz2) < 1e-6, 'Bad derivative');


%% Test 2: Hyperbolic Case

% Choose z > 0 for ellipse
z = -1.5*pi;
dz = 1e-4;
tol = 1e-2;
K = 16;

% Function at z
C = C_func(z, tol, K);
S = S_func(z, tol, K);

% Numerical 1st derivative
C_p = C_func(z+dz, tol, K);
C_m = C_func(z-dz, tol, K);
S_p = S_func(z+dz, tol, K);
S_m = S_func(z-dz, tol, K);
dCdz_true = 0.5*(C_p - C_m)/dz;
dSdz_true = 0.5*(S_p - S_m)/dz;

% Analytical Derivative
dCdz = dCdz_func(z, C, S, tol, K);
dSdz = dSdz_func(z, C, S, tol, K);

% Compare analytical and numerical derivative
assert(abs(dCdz_true - dCdz) < 1e-6, 'Bad derivative');
assert(abs(dSdz_true - dSdz) < 1e-6, 'Bad derivative');

% Numercial 2nd derivative
dCdz_p = dCdz_func(z+dz, C_p, S_p, tol, K);
dCdz_m = dCdz_func(z-dz, C_m, S_m, tol, K);
dSdz_p = dSdz_func(z+dz, C_p, S_p, tol, K);
dSdz_m = dSdz_func(z-dz, C_m, S_m, tol, K);
d2Cdz2_true = 0.5*(dCdz_p - dCdz_m)/dz;
d2Sdz2_true = 0.5*(dSdz_p - dSdz_m)/dz;

% Analytical 2nd derivative
d2Cdz2 = d2Cdz2_func(z, C, S, dCdz, dSdz, tol, K);
d2Sdz2 = d2Sdz2_func(z, C, S, dCdz, dSdz, tol, K);

% Compare analytical and numerical second derivative
assert(abs(d2Cdz2_true - d2Cdz2) < 1e-6, 'Bad derivative');
assert(abs(d2Sdz2_true - d2Sdz2) < 1e-6, 'Bad derivative');


%% Test 3: Parabolic Case at 0

% Choose z = 0 for paraboclic
z = 0;
dz = 1e-4;
tol = 1e-2;
K = 16;

% Function at z
C = C_func(z, tol, K);
S = S_func(z, tol, K);

% Numerical 1st derivative
C_p = C_func(z+dz, tol, K);
C_m = C_func(z-dz, tol, K);
S_p = S_func(z+dz, tol, K);
S_m = S_func(z-dz, tol, K);
dCdz_true = 0.5*(C_p - C_m)/dz;
dSdz_true = 0.5*(S_p - S_m)/dz;

% Analytical Derivative
dCdz = dCdz_func(z, C, S, tol, K);
dSdz = dSdz_func(z, C, S, tol, K);

% Compare analytical and numerical derivative
assert(abs(dCdz_true - dCdz) < 1e-6, 'Bad derivative');
assert(abs(dSdz_true - dSdz) < 1e-6, 'Bad derivative');

% Numercial 2nd derivative
dCdz_p = dCdz_func(z+dz, C_p, S_p, tol, K);
dCdz_m = dCdz_func(z-dz, C_m, S_m, tol, K);
dSdz_p = dSdz_func(z+dz, C_p, S_p, tol, K);
dSdz_m = dSdz_func(z-dz, C_m, S_m, tol, K);
d2Cdz2_true = 0.5*(dCdz_p - dCdz_m)/dz;
d2Sdz2_true = 0.5*(dSdz_p - dSdz_m)/dz;

% Analytical 2nd derivative
d2Cdz2 = d2Cdz2_func(z, C, S, dCdz, dSdz, tol, K);
d2Sdz2 = d2Sdz2_func(z, C, S, dCdz, dSdz, tol, K);

% Compare analytical and numerical second derivative
assert(abs(d2Cdz2_true - d2Cdz2) < 1e-6, 'Bad derivative');
assert(abs(d2Sdz2_true - d2Sdz2) < 1e-6, 'Bad derivative');


%% Test 4: Parabolic Case at z != 0

% Choose z ~= 0 for parabola
z = 2.5e-5;
dz = 1e-4;
tol = 5e-2;
K = 16;

% Function at z
C = C_func(z, tol, K);
S = S_func(z, tol, K);

% Numerical 1st derivative
C_p = C_func(z+dz, tol, K);
C_m = C_func(z-dz, tol, K);
S_p = S_func(z+dz, tol, K);
S_m = S_func(z-dz, tol, K);
dCdz_true = 0.5*(C_p - C_m)/dz;
dSdz_true = 0.5*(S_p - S_m)/dz;

% Analytical Derivative
dCdz = dCdz_func(z, C, S, tol, K);
dSdz = dSdz_func(z, C, S, tol, K);

% Compare analytical and numerical derivative
assert(abs(dCdz_true - dCdz) < 1e-6, 'Bad derivative');
assert(abs(dSdz_true - dSdz) < 1e-6, 'Bad derivative');

% Numercial 2nd derivative
dCdz_p = dCdz_func(z+dz, C_p, S_p, tol, K);
dCdz_m = dCdz_func(z-dz, C_m, S_m, tol, K);
dSdz_p = dSdz_func(z+dz, C_p, S_p, tol, K);
dSdz_m = dSdz_func(z-dz, C_m, S_m, tol, K);
d2Cdz2_true = 0.5*(dCdz_p - dCdz_m)/dz;
d2Sdz2_true = 0.5*(dSdz_p - dSdz_m)/dz;

% Analytical 2nd derivative
d2Cdz2 = d2Cdz2_func(z, C, S, dCdz, dSdz, tol, K);
d2Sdz2 = d2Sdz2_func(z, C, S, dCdz, dSdz, tol, K);

% Compare analytical and numerical second derivative
assert(abs(d2Cdz2_true - d2Cdz2) < 1e-6, 'Bad derivative');
assert(abs(d2Sdz2_true - d2Sdz2) < 1e-6, 'Bad derivative');


%% Ouput
fprintf('PASSED: Test_dCdz\n');