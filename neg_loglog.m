function h = neg_loglog(x,y,h, varargin)

% h = neg_loglog(x,y,h, varargin)
% 
% neg_loglog.m is part of the CR1Dmod forward modeling package, and contains
% the code used to plot data-sets with negative values on a log-log scale.
% For example, the CR-effect in TEM calculations may cause the response in
% certain time ranges to become negative. The negative data are plotted
% as positive, but with a different line type.
%
% Inputs:
% x:         Vector of x values
% y:         Vector of y values
% h:         Uses the figure with handle h for plotting
% varargin:  Additional arguments are passed directly to loglog
%
% Output:
% h:         handle to the plot figure.
%
% Written by:
% Thomas Ingeman-Nielsen
% The Arctic Technology Center, BYG
% Technical University of Denmark
% Email: tin@byg.dtu.dk

if ((isappdata(0,'compiled') && getappdata(0,'compiled')) || exist('h','var')) && ...
        ishandle(h)
    figure(h);
else
    h = figure;
end

if nargin > 3
    plotProp = varargin(1:end);
else
    plotProp = {};
end

if all(isreal(y))
    % in case of real input, make 1 plot
    
    % define segments where all data have the same sign
    segment = 1;
    a = sign(y(1));
    if a==0, a=1; end; 
    for k = 2:length(y)
        a2 = sign(y(k));
        if a2==0, a2=1; end;
        if a2~=a
            segment = [segment k];
        end
        a = a2;
    end
    segment = [segment, length(y)+1];

    % determine the plot properties    
    if sign(y(segment(1)))==1 || sign(y(segment(1)))==0
        color = '-ob-or';
    else
        color = '-or-ob';
    end
    
    % Do the plotting
    for k = 1:length(segment)-1
        loglog(x(segment(k):segment(k+1)-1),y(segment(k):segment(k+1)-1),...
            color(1:3), plotProp{:});
        hold on;
        color = color([4 5 6 1 2 3]);
        y = -y;
    end
else
    % in case of complex input, make 2 plots
    % Treat the real data first
    
    % define segments where all data have the same sign    
    subplot(2,1,1);
    segment = 1;
    a = sign(real(y(1)));
    if a==0, a=1; end;
    for k = 2:length(y)
        a2 = sign(real(y(k)));
        if a2==0, a2=1; end;
        if a2~=a
            segment = [segment k];
        end
        a = a2;
    end
    segment = [segment, length(y)+1];
    
    % determine the plot properties    
    if sign(real(y(segment(1))))==1 || sign(real(y(segment(1))))==0
        color = '-ob-or';
    else
        color = '-or-ob';
    end
    
    % Do the plotting
    for k = 1:length(segment)-1
        loglog(x(segment(k):segment(k+1)-1),abs(real(y(segment(k):segment(k+1)-1)),...
            color(1:3), plotProp{:}));
        hold on;
        color = color([4 5 6 1 2 3]);
    end
    
    % Now treat the imaginary part
    % define segments where all data have the same sign    
    subplot(2,1,2);
    segment = 1;
    a = sign(imag(y(1)));
    if a==0, a=1; end;
    for k = 2:length(y)
        a2 = sign(imag(y(k)));
        if a2==0, a2=1; end;
        if a2~=a
            segment = [segment k];
        end
        a = a2;
    end
    segment = [segment, length(y)];

    % determine the plot properties    
    if sign(imag(y(segment(1))))==1 || sign(imag(y(segment(1))))==0
        color = '-ob-or';
    else
        color = '-or-ob';
    end
    
    % Do the plotting
    for k = 1:length(segment)-1
        loglog(x(segment(k):segment(k+1)-1),abs(imag(y(segment(k):segment(k+1)-1)),...
            color(1:3), plotProp{:}));
        hold on;
        color = color([4 5 6 1 2 3]);
    end
end