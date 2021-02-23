function [sum,i,termstor,eulterm] = NJCST(htype, kernelfcn, r, tol, err, maxite, varargin)

% [sum,i,termstor,eulterm] = NJCST(htype, kernelfn, r, tol, err, maxite, varargin)
%
% NJCST.m is part of the CR1Dmod forward modeling package. It calculates hankel
% or harmonic transform of the the type:
%
%  G(r) = integral of [ FUNC(x) * JCS(x*r) * dx ] from 0 to infinity
%
% where JCS is one of the bessel or harmonic functions J0, J1, COS or SIN.
% The kernel function of the transformation (FUNC) is passed through the 
% function handle 'kernelfcn'. Additional parameters passed to NJCST in the
% varargin cell array, will be passed on to FUNC after the integration 
% parameter x.
% 
% Input
% htype:      'J0', 'J1', 'COS' or 'SIN'
% kernelfcn:  pointer to or expression for kernel function
% r:          necessary constant of the transformation
% maxite:     maximum number of segments of the besselfunction to integrate
%             over (default and absolute max is 100)
% err:        stop integration when last integration result was less than 
%             err.*sum_so_far (for both real and imag part) (default is 1e-5)
% tol:        tolerance for the quadrature routine (default is 1e-6)
%
% Output
% sum:        value of the transformation
% i:          number of segments integrated
% termstor:   array of calculated integration terms
% eulterm:    array of euler-transformed terms. The result of summing all terms
%             in this array, is the output parameter sum
%
% The function uses a quadrature rule to do the integration between pre-calculated 
% zeros of the transform function. For the J0 and J1 functions the first 301
% zeros have been established. If more terms are, these are found by adding
% multiples of pi to the 301st zero. The terms calculated in this manner are
% stored in the termstor array, and an euler transformation is applied. New
% terms are added until the last (euler) term added to the sum is less than
% the sum multiplied by the parameter 'err', or until 'maxite' terms have been 
% added.
% The euler transformation is not applied to the first 'startEuler-1' terms,
% and these terms will be calculated and summed regardless of the error
% criteria.
%
% The Euler transformation code used is modified from:
% Press, W.H., Flannery, B.P., Teukolsky, S.A., and Vetterling, W.T. (1986):
% Numerical Recipes, The Art of Numerical Computing. Cambridge University
% Press.
%
% The calling function can save the last two output terms, termstor and eulterm,
% and use them for calculating related transforms. I.e. it is often necessary
% to calculate another transform where the only difference in the kernel
% function is a power of the integration variable, x. This can be easily
% done by multiplying all the terms in termstor by x and then perform a new
% euler transformation on the terms. However no automatic procedure has been
% implemented to facilitat this, so far.
%
% Written by:
% Thomas Ingeman-Nielsen
% The Arctic Technology Center, BYG
% Technical University of Denmark
% Email: tin@byg.dtu.dk

persistent S  % S will hold precalculated values of the zeros of the transform functions 
startEuler = 3;

% load precalculated zero-values and add more iff necessary.
switch htype
    case {'J0'}
        if ~isfield(S,'z_J') || isempty(S.z_J)
            S = load('Zeros_J.mat','z_J');
        end
        if length(S.z_J)-1<maxite
            S.z_J = [S.z_J, S.z_J(end)+(1:maxite-length(S.z_J)+1).*pi];
        end
        z = S.z_J;
        quadfun = @J0;
    case {'J1'}
        if ~isfield(S,'z_J') || isempty(S.z_J)
            S = load('zeros_J.mat','z_J');
        end
        if length(S.z_J)-1<maxite
            S.z_J = [S.z_J, S.z_J(end)+(1:maxite-length(S.z_J)+1).*pi];
        end
        z = S.z_J;
        quadfun = @J1;
    case {'COS'}
        if ~isfield(S,'z_C') || isempty(S.z_C)
            S.z_C = [0, (1:2:maxite*2+2).*pi./2];
        end
        if length(S.z_C)-1<maxite
            S.z_C = [S.z_C, S.z_C(end)+(1:maxite-length(S.z_C)+1).*pi./2];
        end
        z = S.z_C;
        quadfun = @COSqfun;
    case {'SIN'}
        if ~isfield(S,'z_S') || isempty(S.z_S)
            S.z_S = [0, (2:2:maxite*2+2).*pi./2];
        end
        if length(S.z_S)-1<maxite
            S.z_S = [S.z_S, S.z_S(end)+(1:maxite-length(S.z_S)+1).*pi./2];
        end
        z = S.z_S;
        quadfun = @SINqfun;            
end

if isempty(err)
    err = 1e-5;
end

termstor = [];
eulterm = [];
sum = 0;
i = 1; % i-1 is the number of segments calculated

lastwarn('');
warnState = warning('off', 'MATLAB:divideByZero');
% Since we integrate from 0, kernel functions which contains the term 1/x
% will give a divide by zero warning at the lower limit of integration. This
% is dealt with in quad by adding eps (the smallest number Matlab can
% handle) to the lower limit of integration.
term = quad(quadfun, z(i)/r, z(i+1)/r, tol, 0, kernelfcn, r, varargin{:});
if ~isfinite(term)
    % if quad fails try quadl instead
    term = quadl(quadfun, z(i)/r, z(i+1)/r, tol, 0, kernelfcn, r, varargin{:});
    if ~isfinite(term)
        disp(['NJCST error: neither quad or quadl was able to integrate '...
                'the specified function']);
        disp('Aborting calculation!');
        termstor = [termstor; term];
        return
    end
end
warning(warnState);

termstor = [termstor; term];
sum = sum + term;
i = i+1;

% Calculate a few terms without euler transformation
% and without taking the error criteria into account.
while i-1<startEuler 
    drawnow;  % This command empties the queue

    term = quad(quadfun, z(i)/r, z(i+1)/r, tol, 0, kernelfcn, r, varargin{:});
    if ~isfinite(term)
        % if quad fails try quadl instead
        term = quadl(quadfun, z(i)/r, z(i+1)/r, tol, 0, kernelfcn, r, varargin{:});
        if ~isfinite(term)
            disp(['NJCST error: neither quad or quadl was able to integrate '...
                    'the specified function']);
            disp('Aborting calculation!');
            termstor = [termstor; term];
            return
        end
    end
    termstor = [termstor; term];
    sum = sum + term;
    i = i+1;
end


if (abs(real(term))>=err*abs(real(sum))) || ((abs(imag(sum)) ~= 0) && (abs(imag(term))>=err*abs(imag(sum))))

    term = quadl(quadfun, z(i)/r, z(i+1)/r, tol, 0, kernelfcn, r, varargin{:});
    if ~isfinite(term)
        % if quad fails try quadl instead
        term = quadl(quadfun, z(i)/r, z(i+1)/r, tol, 0, kernelfcn, r, varargin{:});
        if ~isfinite(term)
            disp(['NJCST error: neither quad or quadl was able to integrate '...
                    'the specified function']);
            disp('Aborting calculation!');
            termstor = [termstor; term];
            return
        end
    end
    termstor = [termstor; term];
    i = i+1;

    % The following Euler convergence improvement code is 
    % modified from 'Numerical Recepies in C++'

    wkspc(1) = term;
    nterm = 1;
    sum = sum + 0.5.*term;
    newterm = sum;
    eulterm = zeros(size(termstor)); 

    % Here we begin the main loop
    while i-1<=maxite && ((abs(real(newterm))>=err*abs(real(sum)))...
            || ((abs(imag(sum)) ~= 0) && (abs(imag(newterm))>=err*abs(imag(sum)))))

        drawnow;  % This command empties the queue
        
        term = quadl(quadfun, z(i)/r, z(i+1)/r, tol, 0, kernelfcn, r, varargin{:});
        if ~isfinite(term)
            % if quad fails try quadl instead
            term = quadl(quadfun, z(i)/r, z(i+1)/r, tol, 0, kernelfcn, r, varargin{:});
            if ~isfinite(term)
                disp(['NJCST error: neither quad or quadl was able to integrate '...
                        'the specified function']);
                disp('Aborting calculation!');
                termstor = [termstor; term];
                return
            end
        end
        
        tmp = wkspc(1);
        wkspc(1) = term;
        for k=1:nterm-1
            dum = wkspc(k+1);
            wkspc(k+1) = 0.5.*(wkspc(k)+tmp);
            tmp = dum;
        end
        wkspc(nterm+1) = 0.5.*(wkspc(nterm)+tmp);
        if (abs(wkspc(nterm+1))) <= (abs(wkspc(nterm)))
            nterm = nterm+1;
            newterm = (0.5.*wkspc(nterm));
            sum = sum + newterm;
        else
            newterm = wkspc(nterm+1);
            sum = sum + newterm;
        end
        eulterm = [eulterm;newterm];        
        i = i+1;
    end
end

i = i-1; % This the number of segments included

if isappdata(0,'debug') && getappdata(0,'debug')
    disp(['R: ' sprintf('%6.3f',r) '   N-terms: ' sprintf('%3.0i',i) '   Hank: ' sprintf('%e +i* %e',real(sum), imag(sum))]);
end

% *********************************************************************
% * Kernel functions
% *********************************************************************

function Q = J0(lambda, kernelfcn, r, varargin)
tmp = besselj(0,lambda.*r);
Q = feval(kernelfcn, lambda, varargin{:}).*tmp;

function Q = J1(lambda, kernelfcn, r, varargin)
tmp = besselj(1,lambda.*r);
Q = feval(kernelfcn, lambda, varargin{:}).*tmp;

function Q = COSqfun(lambda, kernelfcn, r, varargin)
tmp = cos(lambda.*r);
Q = feval(kernelfcn, lambda, varargin{:}).*tmp;

function Q = SINqfun(lambda, kernelfcn, r, varargin)
tmp = sin(lambda.*r);
Q = feval(kernelfcn, lambda, varargin{:}).*tmp;




