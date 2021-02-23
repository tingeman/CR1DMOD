function [CCONV, NOUT, ROUT] = FCST(FTYPE, RLO, RHI, IKEEP, IPOW, EPS, FUNC, varargin)

% [CCONV, NOUT, ROUT] = FCST(FTYPE, RLO, RHI, IKEEP, IPOW, EPS, FUNC, varargin);
%
% FCST.m is part of the CR1Dmod forward modeling package. The routine is 
% identical to FJCST, except that it only calculates COS and SIN transforms.
% When calculating transient fields, a Hankel transform occurs in the kernel
% function FUNC, causing problems with changes in the persistent variables, 
% which confuses the SIN/COS transform. The problem doesn't occur when a 
% spline function is used in the frequency domain, since all the frequency 
% domain computations are done before the SIN/COS transform is invoked.
% To avoid the problem, when not using the spline function, two separate
% functions are called for Hankel and Harmonic transforms.
%
% For information about the code and parameters, please refer to FJCST.m
%
% Translated and modified by:
% Thomas Ingeman-Nielsen
% The Arctic Technology Center, BYG
% Technical University of Denmark
% Email: tin@byg.dtu.dk

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
    disp('PARAMETER RLO IN SUBROUTINE FCST OUT OF RANGE.');
    disp([' RLO HAS BEEN SET TO THE MINIMUM VALUE = ' num2str(RLOMIN)])
end

if (RHI > RHIMAX)
    RHI = RHIMAX;
    disp('PARAMETER RHI IN SUBROUTINE FCST OUT OF RANGE.');
    disp([' RHI HAS BEEN SET TO THE MAXIMUM VALUE = ' num2str(RHIMAX)]);
end

if ((IKEEP < 0) || (IKEEP > 1))
    disp('PARAMETER IKEEP IN SUBROUTINE FCST OUT OF RANGE');
end

if ((IKEEP == 1) && (strcmpi(FTYPEOLD, 'XXX')))
    disp('IKEEP=1 AT FIRST CALL TO FCST IS ILLEGAL');
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
    if ismember(FTYPE, {'CO1';'CO2';'CO4';'SI1';'SI2';'SI4'})
        S = load('filters.mat', FTYPE);    % Define filter-coefficients
        FILT(NHLO+201:NHHI+201) = S.(FTYPE)(NHLO+201:NHHI+201);
        FTYPEOLD = FTYPE;
    else
        disp(' CHARACTER STRING "FTYPE" ILLEGAL IN CALL TO FCST');
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

% --- Since COS or SIN filters, multiply with SQRT(PI/2)

for I=1:NOUT
    CCONV(I)=SQPI2*CCONV(I);
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
