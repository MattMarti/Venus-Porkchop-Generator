%% Computes the test for Orbital Elements
% 
% This function is necessary to debug the p-iteration method


%% GPS satellite test

% Data
mu    =   398600.5e9;
toe   =   4.103840000000D+05;
a_true     =   (5.153668163300D+03)^2;
e_true     =   1.446213386950D-01;
M0_true    =   5.830431375610D-01;
omega_true =   1.105691738600D+00;
i0_true    =   9.295239119100D-01;
OMEGA_true =   2.922381628170D+00;

% Generate input
tepoch = toe;
thist = toe + (0:60:3600*12);
ephem_true = [a_true;e_true;M0_true;omega_true;i0_true;OMEGA_true];

% Generate orbit history
[xrvhist_true] = satpos(mu, ephem_true, tepoch, thist);

% Determine orbit parameters
x = xrvhist_true(:,1);
[a, e, M0, omega, i, Omega, nu] = orbital_elements(mu, x);

% Determine cosnu
cosnu = cos(nu);
eccpcosnu = e + cosnu;
onepecccosnu = 1 + e*cosnu;
E0_c = acos(eccpcosnu / onepecccosnu);
M0_c = E0_c - e .* sin(E0_c);

% Determine sinnu
sinnu = sin(nu);
asqrtonemesq = a_true*sqrt(1 - e^2);
r = norm(x(1:3));
E0_s = asin( sinnu * r / asqrtonemesq );
M0_s = E0_s - e .* sin(E0_s);

% Assertion
p = 1e-12;
assert(abs(a - a_true) < a_true*p, 'failed');
assert(abs(e - e_true) < p, 'failed');
assert(abs(omega - omega_true) < p, 'failed');
assert(abs(i - i0_true) < p, 'failed');
assert(abs(Omega - OMEGA_true) < p, 'failed');
assert(abs(M0_true - M0_c) < p, 'failed');
assert(abs(M0_true - M0_s) < p, 'failed');
assert(abs(M0_true - M0) < p, 'failed');

% Regenerate xrvhist and test
ephem = [a;e;M0;omega;i;Omega];
[xrvhist] = satpos(mu, ephem, tepoch, thist);

% Assert
p = 1e-6;
for i = 1:6
    for j = 1:length(thist)
        assert(abs(xrvhist(i,j) - xrvhist_true(i,j)) < p, 'failed')
    end
end

% Out
fprintf('PASSED: Test_orbital_elements\n');

