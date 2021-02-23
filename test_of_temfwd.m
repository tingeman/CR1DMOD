
% EXAMPLE of usage: (copy to command window, without leading "% ")
%
% Calculate the response from a central loop sounding. The transmitter loop
% has a radius of 10 m, the receiver loop an effective area of 100 sqm.
% The ground is two-layer model with the parameters rho_1 = 100 Ohm.m, 
% h1 = 10 m, rho_2 = 30 Ohm.m.
% Layer 1 has a dispersive resistivity of the cole-cole type.
%
% % Setup input structures:
%
config.type = 'TEM Central Loop'  % TEM Central loop configuration
config.TxR  = 7.13                % Transmitter loop radius (m)
%  If using a square transmitter loop, one can specify the side length, TxS,
%  instead, and the function will calculate the radius of the equivalent 
%  circular loop.
%  config.TxS = 40                   % Square transmitter loop side length (m)
config.RxA  = 100                 % Receiver loop effective area (sqm)

layers(1).depth_to_top = 0;       % Layer 1 depth-to-top of layer
layers(1).thickness    = 10;      % Layer 1 thickness
layers(1).rho          = 50;     % Layer 1 resistivity
layers(1).m            = 0.8;     % Layer 1 cole-cole chargeability
layers(1).tau          = 1e-5;    % Layer 1 cole-cole time constant
layers(1).c            = 0.3;     % Layer 1 cole-cole exponent
layers(1).eps_r        = 1;       % Layer 1 relative permittivity
layers(1).mu           = 0;       % Layer 1 magnetic permeability
layers(2).depth_to_top = 10;      % Layer 2 depth-to-top of layer
layers(2).thickness    = Inf;     % Layer 2 thickness (use Inf for bottom layer)
layers(2).rho          = 10;      % Layer 2 resistivity
layers(2).m            = 0;       % Layer 2 cole-cole chargeability
layers(2).tau          = 0;       % Layer 2 cole-cole time constant
layers(2).c            = 0;       % Layer 2 cole-cole exponent
layers(2).eps_r        = 1;       % Layer 2 relative permittivity
layers(2).mu           = 0;       % Layer 2 magnetic permeability

cparams.domain    = 'TD'          % Time Domain calculation
cparams.calc_type = 'Quasi'       % Use quasi-stationary approximation
cparams.times     = logspace(-6,-3,30) % Specify decay times for which to calculate result
cparams.waveform  = 'Step function'    % Only option at the moment! 
cparams.hank_type = 'FHT'         % Choose Fast Hankel Transform using digital filters for calculating FD response
cparams.FHT_err = 1.0000e-08      % Tolerance of Fast Sine Transform
cparams.FDspline = 1              % Use spline to approximate frequency domain response (0=no, 1=yes)
cparams.FDsp_NDEC = 10            % Number of frequencies to calculate for spline interpolation
cparams.FDsp_Bmin = 1.0000e-04    % Smallest frequency (induction number) to calculate
cparams.FDsp_Bmax = 1000          % Largest frequency (induction number) to calculate
cparams.FtoTtype  = 'FST'         % Choose Fast Sine Transform using digital filters
cparams.FST_err   = 1.0000e-08    % Tolerance of Fast Sine Transform
cparams.showFDsp  = 1             % plot underlying freq. domain calculations (0=no, 1=yes)

% The following is only needed for Numerical Hankel Transform integration:

% cparams.Seg_tol = 1e-6;     % Tolerance on each segment of the integration
% cparams.NHT_tol = 1e-5;     % Tolerance on the sum of the series
% cparams.Max_seg = 100;      % Maximum number of segments to sum
% 
% Run forward calculation:

result = temfwd(config,layers,cparams);
plotdat('tem', result, config, cparams)

% END OF EXAMPLE

