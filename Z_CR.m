function Z = Z_CR(omega, R0, tau, c, m, a)

% Z = Z_CR(omega, R0, tau, c, m, a);
%
% Z_CR.m is part of the CR1Dmod forward modeling package, and contains
% the code used to calculate the Cole-Cole resistivity dispersion.
%
% Input parameters:
% omega:       Angular frequency (2.*pi.*freq)
% R0:          Resistivity in the DC limit
% tau:         Time constant
% c:           Frequency dependence
% m:           Chargeability
% a:           Alternative frequency dependence of the 
%                 Cole-Davidson model (optional)
%
% Accepts row-vectors for all inputs (R0, tau, c and m must have 
% same lengths) and returns a matrix with the same number of rows
% as the length of R0, tau, c, and m and the same number of columns 
% as in omega, containing the calculated impedances.
%
% This function may be modified to allow for other types of
% resistivity dispersions.
%
% Written by:
% Thomas Ingeman-Nielsen
% The Arctic Technology Center, BYG
% Technical University of Denmark
% Email: tin@byg.dtu.dk

if nargin < 6
    a = ones(size(c));
end

% Cole-Cole / Cole-Davidson dispersion
Z = R0(:,ones(1,length(omega))).*(1-m(:,ones(1,length(omega)))...
    .*(1-1./(1+(i.*tau*omega).^c(:,ones(1,length(omega))))...
    .^a(:,ones(1,length(omega))))); %impedance 


