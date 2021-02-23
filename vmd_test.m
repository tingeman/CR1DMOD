function vmd_test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vmd.m is included in the CR1Dmod forward modeling package.
% 
% The vmd function illustrates the use of the numerical and fast Hankel
% transform routines of CR1Dmod (FJCST.m and NJCST.m) by calculating the
% response of a Vertical Magnetic Dipole on the surface of a homogeneous
% half-space with a real and frequency independent resistivity.
%
% Written by:
% Thomas Ingeman-Nielsen
% The Arctic Technology Center, BYG
% Technical University of Denmark
% Email: tin@byg.dtu.dk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      
eps = 1./(299792458.^2.*4e-7.*pi);% 8.85418782e-12;   % 1/(c^2.*mu)
mu = 4e-7.*pi;      
disp(' ');
disp('-------------------------------------------------------------------------');
disp('The function calculates the response of a vertical magnetic');
disp('dipole on the surface of a homogeneous half-space.');
disp('The function calculates the response using the fast and numerical');
disp('Hankel transform routines, FJCST.m and NJCST.m, included in CR1Dmod.');
disp('For comparison, the analytical solution is calculated as well.');

% Input half-space resistivity
disp(' ');
rho_in = input('Enter half-space resistivity [default: 100 Ohmm]: ','s');
rho = str2double(rho_in);
while (~isempty(rho_in) && isnan(rho)) || rho<=0
    disp('Bad input, please try again!');
    rho_in = input('Enter half-space resistivity [default: 100 Ohmm]: ','s');
    rho = str2double(rho_in);
end

if isempty(rho_in)
    rho = 100;
end

% Input frequency range
disp(' ');
disp('Enter the frequency or frequency range, using normal Matlab expressions.');
freq_in = input('[default: logspace(2,4,10)]: ','s');
freq = str2num(freq_in);
while (~isempty(freq_in) && isempty(freq)) || any(freq<=0)
    disp('Bad input, please try again!');
    freq_in = input('Enter the frequency or frequency range: ','s');
    freq = str2num(freq_in);
end

if isempty(freq_in)
    freq = logspace(2,4,10);
end

% Input dipole separation
disp(' ');
r_in = input('Enter dipole separation [default: 100 m]: ','s');
r = str2double(r_in);
while (~isempty(r_in) && isnan(r)) || r<=0
    disp('Bad input, please try again!');
    r_in = input('Enter dipole separation [default: 100 m]: ','s');
    r = str2double(r_in);
end

if isempty(r_in)
    r = 100;
end


% Calculate analytical solution

omega = 2.*pi.*freq;
sigma = 1./rho;

k1_sq = omega.^2.*mu.*eps - i.*omega.*mu.*sigma;
k = sqrt(k1_sq);
ana = 1./(2.*pi.*k1_sq.*(r.^5))                          ...
    .*(9 - ( 9 + 9.*i.*sqrt(k1_sq).*r - 4.*k1_sq.*(r.^2)    ...
    - i.*sqrt(k1_sq).^3.*(r.^3) )                        ...
    .*exp(-i.*sqrt(k1_sq).*r) );
   
% To get Z/Z0
free_space = -1./(4.*pi.*(r.^3));
ana = ana./free_space;

% Calculate Fast Hankel Transform

FHT_ans = zeros(size(omega));
tic
for k = 1:length(omega)
    [FHT_ans(k), NOUT, ROUT] = FJCST('J04', r, -1, 0, 1, 1e-8, @vmd_kernel, k1_sq(k));
end
FHT_ans = 1./(2.*pi).*FHT_ans./free_space;
FHT_time = toc;

% Calculate Numerical Hankel Transform

NHT_ans = zeros(size(omega));
nbsg = zeros(size(omega));
tic
for k = 1:length(omega)
    [NHT_ans(k), nbsg(k)] =  NJCST('J0', @vmd_kernel, r, 1e-14, 1e-07, 300, k1_sq(k));
end
NHT_ans = 1./(2.*pi).*NHT_ans./free_space;
NHT_time = toc;

% End of calculations

disp(' ');
disp('Results will be displayed for the first frequency only!');
disp('-------------------------------------------------------');
disp(' ');
disp('                     Real part        Imaginary part');
disp(sprintf('Analytical result:   %.8e  %.8e', real(ana(1)), imag(ana(1))));
disp(sprintf('NHT result:          %.8e  %.8e', real(NHT_ans(1)), imag(NHT_ans(1))));
disp(sprintf('FHT result:          %.8e  %.8e', real(FHT_ans(1)), imag(FHT_ans(1))));
disp(sprintf('NHT nb. of segments: %d', nbsg(1)));
disp(' ');
disp('Calculation times (the sum for all specified frequencies:');
disp(sprintf('NHT time:            %.4d s', NHT_time));
disp(sprintf('FHT time:            %.4d s', FHT_time));
disp(' ');
disp('-------------------------------------------------------------------------');
disp(' ');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction: vmd_kernel
% Kernel function of the Hankel transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kern = vmd_kernel(lambda, k1_sq)
  
u1 = sqrt(lambda.^2 - k1_sq);
kern = (lambda.^3)./(lambda+u1);
