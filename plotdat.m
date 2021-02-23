function plotdat(type, result, config, cparams, varargin)

% Usage:
% plotdat(type, result, config, cparams, fighandle, 'reuse')
%
% plotdat.m is part of the CR1Dmod forward modeling package. The function
% is used to plot the different responses that CR1Dmod can model.
%
% Type can be one of {'emgsa', 'dcgsa', 'hcp', 'tem'}
% if fighandle is not supplied or is not a handle, a new figure is produced
% if the last option is omitted, the figure is overwritten 
%
% Written by:
% Thomas Ingeman-Nielsen
% The Arctic Technology Center, BYG
% Technical University of Denmark
% Email: tin@byg.dtu.dk

if nargin>4 && ishandle(varargin{1})
    fig = varargin{1};
else
    fig = figure('visible', 'off');
end

if nargin == 6 && strcmpi(varargin{2},'reuse')
    children = get(fig,'children');
    ax = findobj(children, 'flat', 'type', 'axes');
    set(ax, 'nextplot', 'add');
end

switch lower(type)
    case{'emgsa'}
        plot_emgsa(result, config, cparams, fig);
    case{'dcgsa'}
        plot_dcgsa(result, config, cparams, fig);
    case{'hcp'}
        plot_hcp(result, config, cparams, fig);
    case{'tem'}
        plot_tem(result, config, cparams, fig);
end

% --------------------------------------------------------------------
function plot_emgsa(result, config, cparams, fig)

mdl.result = result;
mdl.config = config;
mdl.cparams = cparams;

if min(size(result)) == 1
    figure(fig);
    set(fig, 'visible', 'on');
    specs = '+*.xsd^v><ph';
    nsp = 12;
    col = 'brgcmykbrgcmykbrgcmyk';
    cnum = 0;
    for k = 1:length(result)
        if ~isfield(result(k), 'G_factor')        
            result(k).G_factor = 0;
            if isfield(result(k), 'Aspac')
                result(k).G_factor = pi.*result(k).Aspac                ...
                    .*(result(k).Nspac)                                 ...
                    .*(result(k).Nspac+1)                               ...
                    .*(result(k).Nspac+2);
            end
            if isempty(result(k).G_factor) || result(k).G_factor == 0
                result(k).G_factor = 1;
            end % if
        end
        if k>nsp*(cnum+1)
            cnum = cnum +1;
        end
        subplot(2,1,1);
        h = semilogx(cparams.freq, abs(result(k).Z.*result(k).G_factor), ...
            ['-' col(cnum+1)]);
        set(h, 'UserData', mdl);
        if length(result)==1
            set(h, 'marker', 'none');
        else
            set(h, 'marker', specs(k-cnum*nsp), 'markersize', 5);            
        end
        hold on;
        if isfield(result, 'name') && ~isempty(result(1).name)
            title(result(1).name);
        end
        
        subplot(2,1,2);
        h = semilogx(cparams.freq,-angle(result(k).Z.*result(k).G_factor).*1000, ...
            ['-' col(cnum+1)]);
        set(h, 'UserData', mdl);
        if length(result)==1
            set(h, 'marker', 'none');
        else
            set(h, 'marker', specs(k-cnum*nsp), 'markersize', 5);            
        end        
        hold on;
    end % k
    subplot(2,1,1);
    ylabel('\rho_a (\Omega\cdot m)');
    hold off;
    subplot(2,1,2);
    ylabel('neg. phase (mrad)');
    xlabel('Freq. (Hz)');
    hold off;         
else
    % if you want results in a 3D matrix for D-D data
    % Z = reshape([result.Z],length(result(1).Z),  ...
    %    size(result,1),size(result,2));
    disp('I don'' know how to plot these results!');
end % if


% --------------------------------------------------------------------
function plot_dcgsa(result, config, cparams, fig)

if min(size(result)) == 1 && max(size(result)) > 1
    figure(fig);
    set(fig, 'visible', 'on');
    
    yval = squeeze(reshape([result.Z], size(result,1), size(result,2)));
    G_factor = squeeze(reshape([result.G_factor], size(result,1),       ...
        size(result,2)));
    
%    for k = 1:length(yval)
k=1;
        if G_factor(k) == 0
            G_factor(k) = 1;
        end % if
        switch config.type
            case {'Dipole-Dipole', '*Capacitance*'}
                if size(result,1) > 1
                    xval = squeeze(reshape(         ...
                        [result.Aspac],             ...
                        size(result,1), 1));
                    loglog(xval, yval.*G_factor);
                else
                    xval = squeeze(reshape(         ...
                        [result.Nspac],             ...
                        1, size(result,2)));
                    semilogy(xval, yval.*G_factor);
                end % if
            case {'Wenner'}
                if size(result,1) > 1
                    xval = squeeze(reshape(             ...
                        [result.Aspac],                 ...
                        size(result,1), 1));
                    semilogy(xval, yval.*G_factor);
                else
                    xval = squeeze(reshape(             ...
                        [result.Aspac],                 ...
                        1, size(result,2)));
                    semilogy(xval, yval.*G_factor);
                end
            case {'Schlumberger'}
                if size(result,1) > 1
                    xval = squeeze(reshape(         ...
                        [result.OA],                ...
                        size(result,1), 1));
                    loglog(xval, yval.*G_factor);                                                                            
                else
                    xval = squeeze(reshape(         ...
                        [result.OM],                ...
                        1, size(result,2)));
                    semilogy(xval, yval.*G_factor);
                end % if
        end
        hold on;
        %end % k
    hold off;
    ylabel('Apparent Resistivity (\Omegam)');
    switch config.type
        case {'Dipole-Dipole', '*Capacitance*'}
            if size(result,1) > 1
                xlabel('A-spacing (m)');
            else
                xlabel('N-spacing (m)');
            end % if
        case {'Wenner'}
            xlabel('A-spacing (m)');
        case {'Schlumberger'}
            if size(result,1) > 1
                xlabel('OA (m)');                                                                            
            else
                xlabel('OM (m)');
            end % if
    end % switch   
end % if

% --------------------------------------------------------------------
function plot_hcp(result, config, cparams, fig)

if min(size(result)) == 1
    fig = figure(fig); 
    set(fig, 'visible', 'on');
    
    for k = 1:length(result)
        subplot(2,1,1);
        semilogx(cparams.freq, real(result(k).H./result(k).H_prim));
        hold on;
        subplot(2,1,2);
        semilogx(cparams.freq, imag(result(k).H./result(k).H_prim));
        hold on;
    end % k
    subplot(2,1,1);
    ylabel('Real(Z/Z_0)');
    hold off;
    subplot(2,1,2);
    ylabel('Imag(Z/Z_0)');
    xlabel('Freq. (Hz)');
    hold off;         
else
    disp('I don'' know how to plot these results!');
end % if


% --------------------------------------------------------------------
function plot_tem(v, config, cparams, fig)

if getappdata(0, 'debug')
    disp('comput2.m:Plot_TEM: Debug before plot...');
end

fig = figure(fig); 
set(fig, 'visible', 'on');
if strcmp(cparams.domain,'TD')
    neg_loglog(cparams.times,v,fig,'Markersize',5);
    ylabel('V/A');
    % rhoa = (config.TxR.^(4/3).*config.RxA.^(2/3).*(4e-7*pi).^(5/3))./ ...
    %     (20.^(2/3).*pi.^(1/3).*cparams.times.^(5/3).*v.^(2/3));
    % loglog(cparams.times,rhoa,'-xb','Markersize',5);
    % ylabel('\rho_a (\Omega m)');
    xlabel('Time (sec)');
else
    subplot(2,1,1);
    semilogx(cparams.freq, real(v), '-xb', 'Markersize',5);
    ylabel('Real(H_z)');
    subplot(2,1,2);
    semilogx(cparams.freq, imag(v), '-xb', 'Markersize',5);
    ylabel('Imag(H_z)');
    xlabel('Frequency (Hz)');
end
