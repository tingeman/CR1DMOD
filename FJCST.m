function [CCONV, NOUT, ROUT] = FJCST(FTYPE, RLO, RHI, IKEEP, IPOW, EPS, FUNC, varargin)

% [CCONV, NOUT, ROUT] = FJCST(FTYPE, RLO, RHI, IKEEP, IPOW, EPS, FUNC, varargin);
%
% FCST.m is part of the CR1Dmod forward modeling package. It is a function used
% to calculate a J0, J1, COS, or SIN transform of a complex or real function.
% The integral is computed as a convolution between sampled values of the 
% kernel function and Fast Hankel Transform filter coefficients.
%
%  What is calculated is the integral:
%  
%  G(r)=integral of [ FUNC(x) * (x**IPOW) * JCS(x*r) * dx ]
%  
%  from zero to infinity, where JCS is either J0, J1, COS, or SIN.
% 
% Input parameters:
% FTYPE:    Type of transform, use 'J04' for electromagnetic problems
% RLO:      Smallest R-value to be calculated
% RHI:      Largest R-value to be calculated, when negative, only RLO
%             is calculated
% IKEEP:    Keep previously calculated kernel function evaluations if 1
% IPOW:     Power of lambda multiplied with kernel function, used only 
%             when IKEEP is 1
% EPS:      Desired relative accuracy
% FUNC:     Function handle to kernel function
% varargin: Cell-array of variables passed directly to FUNC
%
% Output parameters:
% CCONV:    Result of the transform, array of values for chosen R-values
% NOUT:     Number of R-values calculated
% ROUT:     The actual R-values used in the calculations
%
% The kernel function must accept an array of values of the transform 
% parameter x as the first input, and any additional inputs should be
% transferred using the 'varargin' cell-array. The output must be an 
% array of values of identical size.
%
% This function has been translated from the cconvol Fortran routine 
% originally written by: 
% Niels Bøie Christensen, 
% Department of Earth Sciences,
% University of Aarhus
%
% For more information about the choice of parameters please refer to 
% the notes from the original routine included below.
%
% Translated and modified by:
% Thomas Ingeman-Nielsen
% The Arctic Technology Center, BYG
% Technical University of Denmark
% Email: tin@byg.dtu.dk
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%     Original comments for CCONVOL                                   %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%     S U B R O U T I N E   C C O N V O L                             %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Subroutine CCONVOL calculates a J0, J1, COS, or SIN transform of a
%  COMPLEX function. The integral is computed as a convolution between
%  sampled values of the kernel function and Fast Hankel Transform
%  filter coefficients. All calculations are in double precision
%  REAL*8 and COMPLEX*16.
%
%  What is calculated is the integral:
%  
%  G(r)=integral of [ FUNC(x) * (x**IPOW) * JCS(x*r) * dx ]
%  
%  from zero to infinity, where JCS is either J0, J1, COS, or SIN.
%  
%  INPUT PARAMETERS:
%  =================
%   FTYPE (CHARACTER*3) determines the type of transform. FTYPE must be
%   given with capital letters.
%   FTYPE='J01' :  J0-transform with omega=pi filters
%                  (analytical input functions)
%   FTYPE='J02' :  J0-transform with omega=pi/2 filters
%                  (e.g. geoelectrical input functions)
%   FTYPE='J04' :  J0-transform with omega=pi/4 filters
%                  (e.g. electromagnetic input functions)
%   J1-transforms are given by FTYPE=J11, J12, and J14, respectively.
%   COS-transforms are given by FTYPE=CO1, CO2, and CO4, respectively.
%   SIN-transforms are given by FTYPE=SI1, SI2, and SI4, respectively.
%
%  The integral is calculated for logarithmically regularly spaced
%  values of R.
%
%  R = exp(N*DEL)  ,  min(N) = NGLO  ,   max(N) = NGHI
%
%  Calculation is performed for values of R between R1 and R2, where R1
%  is the second largest regular R-value smaller than RLO, and R2 is
%  the second smallest regular R-value greater than RHI.
%
%  If RHI is negative computation is done for R=RLO only and RLO does
%  not have to be one of the regular R-values above, but can have any
%  value. This option makes it possible to avoid interpolation, if only
%  one R-value is needed, e.g. frequency soundings with only one
%  transmitter/receiver separation.
%
%  FUNC is the COMPLEX kernel function.
%
%  IKEEP (INTEGER) determines if the kernel function values are to be kept from
%  the previous computation.
%  IKEEP=1 :  keep the old values
%  IKEEP=0 :  compute new values
%  If IKEEP=1 then the kernel function FUNC must be the same as in the
%  previous call, but the input function FUNC is now multiplied with
%  a power of the integration variable.
%  
%  IPOW (INTEGER) the power of the integration variable with which the
%  input function is multiplied, when IKEEP=1. Only the values
%  IPOW=-1, 1, 2 are allowed. If other powers or other functions are
%  needed, please use the convolution routine with the option of
%  including any user defined function.
%  This option has been implemented to save computation time for
%  related transforms often met in e.g. EM calculations.
%  
%  EPS is the desired relative accuracy of the calculation.
%
%  If RLO, or RHI violate their restricted intervals determined by the
%  dimensioning of the subroutine, they are reset to their limiting
%  values.
%  
%  OUTPUT PARAMETERS:
%  ==================
%  NOUT is the number of values of ROUT, where the integral is
%  calculated.
%  
%  CCONV is a COMPLEX array containing the values of the integral.
%  
%  ROUT and CCONV must be dimensioned:
%  REAL*8 ROUT(NG)
%  COMPLEX*16 CCONV(NG)
%  in the calling subroutine, where NG must be the same as in this
%  subroutine.
%  
%  Subroutine CCONVOL needs a block data subprogram to transfer the
%  filter coefficients for the Hankel transform through the common
%  blocks /J0BLOCK/ and /J1BLOCK/.
%  
%  DIMENSIONING OF THE ARRAYS:
%  ===========================
%  The subroutine is dimensioned through the PARAMETER statements in
%  the first few lines of the code.
%
%  The filter coefficients used in the subroutine are from number NHLO 
%  to number NHHI.
%  
%  The maximum interval within which the output function can be
%  calculated is given by the numbers NGLO and NGHI (see above under
%  input parameters).
%  
%  The sampling density of the filters used is given by the statement
%  "NDEC=??" indicating the number of samples per decade. The sampling
%  density can also be given by the cut-off frequency (wavenumber) SC,
%  in which case the above statement should be changed to "SC=??".
%  
%  NLIM is the number of terms in the discrete convolution, which are
%  calculated initially without checking if the desired relative accuracy
%  has been reached. Recommended value is NLIM=4*NDEC.
%  This saves computation time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent FTYPEOLD KMIN KMAX FILT FC
       
NHLO = -200;
NHHI = 100;
NGLO= -20;
NGHI = 50;
%NG = NGHI-NGLO+1;
NFLO = NGLO-NHHI;
NFHI = NGHI-NHLO;
NDEC = 10;
NLIM = 4*NDEC;

% ----------

if isempty(FTYPEOLD)
    FTYPEOLD = 'XXX';
    FILT = zeros(NHHI-NHLO+1,1);  % will store filter coefficients
%    FC = zeros(NFHI-NFLO+1,1);    % This line moved to test of KEEP condition    
end

SQPI2 = 1.2533141373155;    % sqrt(pi/2)

% ----------

KEEP = IKEEP;
IPOT = IPOW;

% --- CALCULATE MAX AND MIN RADIAL DISTANCE FROM INITIAL PARAMETERS

if NDEC ~= 0
    SC = NDEC/log(100);
else
    disp('Error: 0 samples per decade defined for NDEC');
    return;
end

DEL = 0.5/SC;
RLOMIN=exp(NGLO*DEL);
RHIMAX=exp(NGHI*DEL);

% --- CHECK FOR ERRONEUS PARAMETERS IN THE SUBROUTINE CALL.

if (RLO < RLOMIN) 
    RLO=RLOMIN;
    disp('PARAMETER RLO IN SUBROUTINE FJCST OUT OF RANGE.');
    disp([' RLO HAS BEEN SET TO THE MINIMUM VALUE = ' num2str(RLOMIN)])
end

if (RHI > RHIMAX)
    RHI = RHIMAX;
    disp('PARAMETER RHI IN SUBROUTINE FJCST OUT OF RANGE.');
    disp([' RHI HAS BEEN SET TO THE MAXIMUM VALUE = ' num2str(RHIMAX)]);
end

if ((IKEEP < 0) || (IKEEP > 1))
    disp('PARAMETER IKEEP IN SUBROUTINE FJCST OUT OF RANGE');
end

if ((IKEEP == 1) && (strcmp(FTYPEOLD, 'XXX')))
    disp('IKEEP=1 AT FIRST CALL TO FJCST IS ILLEGAL');
    disp(' PROGRAM STOPPED!');
    return
end    

if ~exist('filters.mat', 'file')
    if ~exist('define_filters.m','file')
        disp(' filters.mat not found, unable to load filter coefficients!');
        disp(' PROGRAM STOPPED!');
        return
    end
    define_filters;
end

% --- THE CHOSEN FILTER COEFFICIENTS ARE PUT INTO THE ARRAY FILT

if ~strcmp(FTYPE,FTYPEOLD)  
    if ismember(FTYPE, {'J01';'J11';'J02';'J12';'J04';...
                'J14';'CO1';'CO2';'CO4';'SI1';'SI2';'SI4'})
        S = load('filters.mat', FTYPE);    % Define filter-coefficients
        FILT(NHLO+201:NHHI+201) = S.(FTYPE)(NHLO+201:NHHI+201);
        FTYPEOLD = FTYPE;
    else
        disp(' CHARACTER STRING "FTYPE" ILLEGAL IN CALL TO FJCST');
        disp(' PROGRAM STOPPED!');
        return
    end
end

% --- INITIALIZATIONS

E = exp(DEL);
E1 = 1/E;

if (RHI < 0)    % If only one R-distance required
    R1 = RLO;
    R2 = RLO;
    NLO = floor(log(RLO)/DEL+100)-100;
    R10 = exp(NLO*DEL);
    XFAC = R10/R1;
    NOUT = 1;
else            % Otherwise...
    NLO = floor(log(RLO)/DEL+100)-101;
    if (NLO < NGLO) 
        NLO = NGLO;
    end
    NHI = floor(log(RHI)/DEL+100)-98;
    if (NHI > NGHI) 
        NHI = NGHI;
    end
    NOUT = NHI-NLO+1;
    R1 = exp(NLO*DEL);
    R2 = exp(NHI*DEL);
    XFAC = 1;
end

% --- Define initial interval if FUNC is new

if (IKEEP == 0)
    KMIN = NLO;
    KMAX = NLO+NLIM;
    
    X = E/R1;
    X = X.*E1.^(1:NLIM+1);

    FC = zeros(NFHI-NFLO+1,1);
    FC(KMIN+201:KMAX+201) = feval(FUNC,X,varargin{:});
end

% --- CALCULATIONS BEGIN HERE

% --- Calculate for smallest R-value

R = R1;

XFIRST = XFAC*E1^(KMIN-1);
S=0;
S = CONVOF(NLO,KMIN,KMAX,1,S,XFIRST,E1,FC,FILT,KEEP,IPOT);

if (S ~= 0)
    XFIRST = XFAC*E1^KMAX;
    [S,FC,KLAST] = CONVON(NLO,KMAX+1,NLO-NHLO,1,S,XFIRST,E1,FUNC,EPS,FC,FILT,KEEP,IPOT,varargin{:});
    KMAX = KLAST;
    
    XFIRST = XFAC*E1^KMIN;
    [S,FC,KLAST] = CONVON(NLO,KMIN-1,NLO-NHHI,-1,S,XFIRST,E,FUNC,EPS,FC,FILT,KEEP,IPOT,varargin{:});
    KMIN = KLAST;
end
CCONV(1) = S/R;
ROUT(1) = R;

if (NOUT ~= 1)
    
    % --- Calculate for greatest R-value
    
    R = R2;
    
    XFIRST = E1^(KMIN-1);
    S=0;
    S = CONVOF(NHI,KMIN,KMAX,1,S,XFIRST,E1,FC,FILT,KEEP,IPOT);
    if (S ~= 0)
        XFIRST = E1^KMAX;
        [S,FC,KLAST] = CONVON(NHI,KMAX+1,NHI-NHLO,1,S,XFIRST,E1,FUNC,EPS,FC,FILT,KEEP,IPOT,varargin{:});
        KMAX = KLAST;
    end
    
    CCONV(NOUT) = S/R;
    ROUT(NOUT) = R;
    
    if (NOUT > 2)
        
        % --- Calculate for all other R-values
        
        R = R1;
        for I=NLO+1:NHI-1
            S=0;
            XFIRST = E1^(KMIN-1);
            S = CONVOF(I,KMIN,KMAX,1,S,XFIRST,E1,FC,FILT,KEEP,IPOT);
            R = R*E;
            CCONV(I-NLO+1) = S/R;
            ROUT(I-NLO+1) = R;
        end
    end
end

% --- If COS or SIN filters, multiply with SQRT(PI/2)


if ~strcmp(FTYPE(1),'J')
    
    for I=1:NOUT
        CCONV(I)=SQPI2*CCONV(I);
    end
    
end

% --- END CALCULATIONS


%CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
%C                                                                     C
%C     S U B R O U T I N E   C O N V O F                               C
%C                                                                     C
%CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function S = CONVOF(IR,K1,K2,KDEL,S,X,DX,FC,FILT,KEEP,IPOT)

%NHLO = -200;
%NHHI = 100;
%NGLO = -20;
%NGHI = 50;
%NG = NGHI-NGLO+1;

%NFLO = NGLO-NHHI;
%NFHI = NGHI-NHLO;

% ----------

K = (K1:KDEL:K2);

if (KEEP == 0)
  
    S=S+sum(FC(K+201).*FILT(IR-K+201));
    
else
    X = X.*DX.^(1:length(K));    
    if (IPOT == -1)
        S = S+sum(FC(K+201).*FILT(IR-K+201)./X);
    end
    
    if (IPOT == 1)
        S = S+sum(FC(K+201).*FILT(IR-K+201).*X);
    end
    
    if (IPOT == 2)
        S = S+sum(FC(K+201).*FILT(IR-K+201).*X.^2);
    end
end


%CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
%C                                                                     C
%C     S U B R O U T I N E   C O N V O N                               C
%C                                                                     C
%CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
function [S,FC,KLAST] = CONVON(IR,K1,K2,KDEL,S,X,DX,FUNC,EPS1,FC,FILT,KEEP,IPOT,varargin)

%NHLO = -200;
%NHHI = 100;
%NGLO = -20;
%NGHI = 50;
%NG = NGHI-NGLO+1;

%NFLO = NGLO-NHHI;
%NFHI = NGHI-NHLO;

% ----------

if (KEEP == 0)
    
    for K=K1:KDEL:K2
        X = X*DX;
        FC(K+201) = feval(FUNC,X,varargin{:});
        SDEL = FC(K+201)*FILT(IR-K+201);
        S = S+SDEL;
        if (abs(SDEL/S) < EPS1)
            break
        end
    end
    
else
    
    if (IPOT == -1)
        for K=K1:KDEL:K2
            X = X*DX;
            FC(K+201) = feval(FUNC,X,varargin{:});
            SDEL = FC(K+201)*FILT(IR-K+201)/X;
            S = S+SDEL;
            if (abs(SDEL/S) < EPS1)
                break
            end
        end
        
    elseif (IPOT == 1)
        for K=K1:KDEL:K2
            X = X*DX;
            FC(K+201) = feval(FUNC,X,varargin{:});
            SDEL = FC(K+201)*FILT(IR-K+201)*X;
            S = S+SDEL;
            if (abs(SDEL/S) < EPS1)
                break
            end
        end
    elseif (IPOT == 2)
        for K=K1:KDEL:K2
            X = X*DX;
            FC(K+201) = feval(FUNC,X,varargin{:});
            SDEL = FC(K+201)*FILT(IR-K+201)*X*X;
            S = S+SDEL;
            if (abs(SDEL/S) < EPS1)
                break
            end
        end
    end
    
end

KLAST = K;
