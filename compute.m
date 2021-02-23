function varargout = compute(varargin)

% compute.m is part of the CR1Dmod forward modeling package.
% 
% It contains functions relating to the compute GUI, from which calculation
% speciffic parameters are defined and the actual computation started. The 
% functions for storing and plotting results are also located in this file.
%
% Default values for the different parameters are located in the subfunction
% "initialize_gui", and may be modified by the user.
%
% Functions for plotting and storing the results are located at the end of
% the file.
%
% Written by:
% Thomas Ingeman-Nielsen
% The Arctic Technology Center, BYG
% Technical University of Denmark
% Email: tin@byg.dtu.dk

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @compute_OpeningFcn, ...
                   'gui_OutputFcn',  @compute_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before compute is made visible.
function compute_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to compute (see VARARGIN)


% Choose default command line output for compute
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if ~isempty(varargin)
    if ishandle(varargin{1})
        handles.modelwin = varargin{1};  % handle to model window
    else
        return;          % bad input
    end
else
    return;              % bad input
end

guidata(hObject, handles);


if strcmp(get(hObject,'Visible'),'off')
    initialize_gui(hObject, handles);
end

handles = guidata(hObject);

setup_gui(hObject, [], handles, 'Default');
% UIWAIT makes compute wait for user response (see UIRESUME)
%uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = compute_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

if isstruct(handles) & isfield(handles,'output')
    varargout{1} = handles.output;
else
    varargout{1} = [];
end

% -------------------------------------------------------------------------
function calculate_batch(hObject, eventdata, handles)

storeAllSettings(hObject, [], handles);

modelwin_han = guidata(handles.modelwin);

model.cparams = modelwin_han.cparams;
model.config = modelwin_han.config;
model.layers = modelwin_han.layers;

batmodels = modelwin_han.batch;

for k = 1:length(batmodels)
    Calculate_button_Callback(hObject, eventdata, handles, batmodels(k));
end
    
% -------------------------------------------------------------------------
function Calculate_button_Callback(hObject, eventdata, handles, varargin)

if nargin == 4 && isstruct(varargin{1})
    model = varargin{1};
    batch = 1;
else
    batch = 0;
    storeAllSettings(hObject, [], handles);
    
    modelwin_han = guidata(handles.modelwin);
    
    model.cparams = modelwin_han.cparams;
    model.config = rmfield(modelwin_han.config,{'plot_handle'});
    model.layers = rmfield(modelwin_han.layers,{'label_handle'});
    if isfield(model.layers,'interface_handle')
        model.layers = rmfield(model.layers,{'interface_handle',            ...
                'prop_lab_handle'});
    end
end

if ~isfield(model,'name')
    model.name = [];
end

for k = 1:length(model.layers)
    if isstr(model.layers(k).mu)
        model.layers(k).mu = str2double(model.layers(k).mu);
    end
end

if model.cparams.domain ~= 'DC'
    if isfield(model.cparams,'times')
        model.cparams.times = str2num(model.cparams.times);
    else
        model.cparams.freq = str2num(model.cparams.freq);
    end
end
result = [];

% prepare plotting
if ~isfield(handles,'fig') || isempty(handles.fig)       || ...
        ~get(handles.Reusefig_checkbox,'value')          || ...
        ~strcmp(model.config.type,handles.figtype{1})    || ...
        ~strcmp(model.cparams.domain,handles.figtype{2}) || ...
        ~ishandle(handles.fig)
    handles.fig = figure('visible','off');
    handles.figtype = {model.config.type;model.cparams.domain};
    guidata(hObject,handles);
else
    children = get(handles.fig,'children');
    ax = findobj(children, 'flat', 'type', 'axes');
    set(ax, 'nextplot', 'add');
    lh = findobj(cell2mat(get(ax,'children')), 'flat', 'type', 'line');
    set(lh, 'marker', 'none');
end

switch model.config.type
    case {'TEM Central Loop'}
        result = temfwd(model.config, model.layers, model.cparams);
        if ~isempty(result)
            %plot_tem(result, model.config, model.cparams, handles.fig);
            plotdat('tem',result, model.config, model.cparams, handles.fig);
        end
    case {'HCP FDEM (HLEM)'}
        result = fdemfwd(model.config, model.layers, model.cparams);
        if ~isempty(result)
            %plot_hcp(result, model.config, model.cparams, handles.fig);
            plotdat('hcp',result, model.config, model.cparams, handles.fig);
        end % if
    otherwise
        switch model.cparams.domain
            case {'FD'}
                switch model.config.type
                    case {'Dipole-Dipole',                             ...
                                'General Surface Array'}
                        result = emgsafwd(model.config,                ...
                            model.layers, model.cparams, model.name);
%                     case {'*Capacitance*', '*GSA capacitance*'}
%                         result = emgsa_cap_fwd(model.config,                ...
%                             model.layers, model.cparams);
                    otherwise
                        disp(['The configuration is not implemented'    ...
                                ' for frequency domain calculations.']);
                        return
                end
                if ~isempty(result)
                    if isfield(model, 'name')
                        [result.name] = deal(model.name);
                    end
                    %plot_emgsa(result, model.config,                    ...
                    %    model.cparams, handles.fig);    
                    plotdat('emgsa', result, model.config,               ...
                        model.cparams, handles.fig);    
                end % if
            case {'DC'}
                result = dcgsafwd(model.config, model.layers,          ...
                    model.cparams);
                if ~isempty(result)
                    plotdat('dcgsa', result, model.config, model.cparams, handles.fig);
                end % if
            otherwise
                disp(['Time domain calculations not yet implemented '   ...
                        'for this configuration']);        
        end % switch
end % switch

if ~isempty(result)
    if getappdata(0, 'debug')
        disp('compute.m:Calculate_button_Callback:debug:');
        disp('    Result variable exists and is not empty!');
    end
    
    if isfield(model.config,'plot_handle')
        rmfield(model.config,'plot_handle');
    end
    for k = 1:length(model.layers)
        if isfield(model.layers(k),'label_handle')
            rmfield(model.layers(k),'label_handle');
        end
        if isfield(model.layers(k),'interface_handle')
            rmfield(model.layers(k),'interface_handle');
        end
        if isfield(model.layers(k),'thick_lab_handle')
            rmfield(model.layers(k),'thick_lab_handle');
        end
    end
    model.result = result;
    % always save to 'forward.mat'
    save('forward.mat','model','-mat');

    % update plots before asking to save
    drawnow;
    if batch
        if ~isfield(model, 'name') || isempty(model.name)
            model.name = 'Model';
        end
        filename = [model.name '_' strrep(datestr(now),':','-')];
        filename = strrep(filename,' ','_');
        saveforward(model, filename);
    else
        saveforward(model);
    end
    %    For standalone mode
    %    assignin('base', 'model',model);
end


% -------------------------------------------------------------------------
function initialize_gui(hObject, handles)

handles.DefaultStrings = {                          ;                   ...
        handles.Times_edit,     'logspace(-6,-3,30)';                   ...
        handles.NST_tol_edit,   '1e-8'              ;                   ...
        handles.FST_err_edit,   '1e-8'              ;                   ...
        handles.Freq_edit,      'logspace(-1,4,25)' ;                   ...
        handles.NHT_tol_edit,   '1e-8'              ;                   ...
        handles.FHT_err_edit,   '1e-8'              ;                   ...
        handles.FDsp_NDEC_edit, '10'                ;                   ...
        handles.FDsp_Bmin_edit, '0.001'             ;                   ...
        handles.FDsp_Bmax_edit, '1000'              ;                   ...
        handles.Quad_tol_edit,  '1e-5'              ;                   ...
        handles.Rsp_NDEC_edit,  '10'                ;                   ...
        handles.Seg_tol_edit,   '1e-12'             ;                   ...
        handles.Max_seg_edit,   '300'               };

for k = 1:length(handles.DefaultStrings);
    set(cell2mat(handles.DefaultStrings(k,1)), 'String',                ...
        handles.DefaultStrings{k,2});
end % k

%PropName   =                      {'Enable', 'Value'};
Settings(1).Config  = {'TEM Central Loop'};
Settings(1).defaultDomain = 3; % = Time Domain
Settings(1).Defaults = {...        %  DC         FD	        TD     group
        handles.DC_radiobutton,      'off', [], 'off', 0 , 'off', 0 , 1 ;...  1
        handles.Freqdom_radiobutton, 'off', [], 'on' , 1 , 'on' , 0 , 1 ;...  2
        handles.Timedom_radiobutton, 'off', [], 'on' , 0 , 'on' , 1 , 1 ;...  3
        handles.Full_sol_radiobutton,'off', [], 'on' , 0 , 'on' , 0 , 2 ;...  4
        handles.Quasi_radiobutton,   'off', [], 'on' , 1 , 'on' , 1 , 2 ;...  5
        handles.Times_edit,          'off', [], 'off', [], 'on' , [], 3 ;...  6
        handles.Times_txt,           'off', [], 'off', [], 'on' , [], 3 ;...  7
        handles.NST_radiobutton,     'off', [], 'off', 0 , 'on' , 0 , 4 ;...  8
        handles.NST_tol_edit,        'off', [], 'off', [], 'off', [], 4 ;...  9
        handles.NST_tol_txt,         'off', [], 'off', [], 'off', [], 4 ;... 10
        handles.FST_radiobutton      'off', [], 'off', 0 , 'on' , 1 , 4 ;... 11
        handles.FST_err_edit,        'off', [], 'off', [], 'on' , [], 4 ;... 12
        handles.FST_err_txt,         'off', [], 'off', [], 'on' , [], 4 ;... 13
        handles.showFDsp_checkbox,   'off', [], 'off', 0 , 'on' , 0 , 5 ;... 14
        handles.Waveform_popup,      'off', [], 'off', 1 , 'off', 1 , 6 ;... 15
        handles.Waveform_txt,        'off', [], 'off', [], 'off', [], 6 ;... 16
        handles.Freq_edit,           'off', [], 'on' , [], 'off', [], 7 ;... 17
        handles.Freq_txt,            'off', [], 'on' , [], 'off', [], 7 ;... 18
        handles.NHT_radiobutton,     'off', [], 'on' , 0 , 'on' , 0 , 8 ;... 19
        handles.NHT_tol_edit,        'off', [], 'off', [], 'off', [], 8 ;... 20
        handles.NHT_tol_txt,         'off', [], 'off', [], 'off', [], 8 ;... 21
        handles.FHT_radiobutton,     'off', [], 'on' , 1 , 'on' , 1 , 8 ;... 22
        handles.FHT_err_edit,        'off', [], 'on' , [], 'on' , [], 8 ;... 23
        handles.FHT_err_txt,         'off', [], 'on' , [], 'on' , [], 8 ;... 24
        handles.FDspline_checkbox,   'off', [], 'off', 0 , 'on' , 1 , 9 ;... 25
        handles.FDsp_NDEC_edit,      'off', [], 'off', [], 'on' , [], 9 ;... 26
        handles.FDsp_NDEC_txt,       'off', [], 'off', [], 'on' , [], 9 ;... 27
        handles.FDsp_Bmin_Bmax_txt,  'off', [], 'off', [], 'on' , [], 9 ;... 28
        handles.FDsp_Bmin_edit,      'off', [], 'off', [], 'on' , [], 9 ;... 29
        handles.FDsp_Bmax_edit,      'off', [], 'off', [], 'on' , [], 9 ;... 30
        handles.Quad_tol_edit,       'off', [], 'off', [], 'off', [], 10;... 31
        handles.Quad_tol_txt,        'off', [], 'off', [], 'off', [], 10;... 32
        handles.Rspline_checkbox,    'off', [], 'off', 0 , 'off', 0 , 11;... 33
        handles.Rsp_NDEC_txt,        'off', [], 'off', [], 'off', [], 11;... 34
        handles.Rsp_NDEC_edit,       'off', [], 'off', [], 'off', [], 11;... 35
        handles.Seg_tol_edit,        'off', [], 'off', [], 'off', [], 12;... 36
        handles.Seg_tol_txt,         'off', [], 'off', [], 'off', [], 12;... 37
        handles.Max_seg_edit,        'off', [], 'off', [], 'off', [], 12;... 38
        handles.Max_seg_txt,         'off', [], 'off', [], 'off', [], 12};%  39   
Settings(1).includeInComparison = find(~ismember(                       ...
    [1:size(Settings(1).Defaults,1)], [9 10 12 13 20 21 23 24 36:39]));

Settings(2).Config  = {'Dipole-Dipole';                                 ...
                       '*Capacitance*';                                 ...
                       'General Surface Array';                         ...
                       '*GSA capacitance*'};
Settings(2).defaultDomain = 1; % = DC
Settings(2).Defaults = {...        %  DC	 FD	    TD	      group
        handles.DC_radiobutton,      'on' , 1 , 'on' , 0 , 'off', [], 1 ;...
        handles.Freqdom_radiobutton, 'on' , 0 , 'on' , 1 , 'off', [], 1 ;...
        handles.Timedom_radiobutton, 'off', 0 , 'off', 0 , 'off', [], 1 ;...
        handles.Full_sol_radiobutton,'off', 0 , 'on' , 0 , 'off', [], 2 ;...
        handles.Quasi_radiobutton,   'off', 0 , 'on' , 1 , 'off', [], 2 ;...
        handles.Times_edit,          'off', [], 'off', [], 'off', [], 3 ;...
        handles.Times_txt,           'off', [], 'off', [], 'off', [], 3 ;...
        handles.NST_radiobutton,     'off', 0 , 'off', 0 , 'off', [], 4 ;...
        handles.NST_tol_edit,        'off', [], 'off', [], 'off', [], 4 ;...
        handles.NST_tol_txt,         'off', [], 'off', [], 'off', [], 4 ;...
        handles.FST_radiobutton      'off', 0 , 'off', 0 , 'off', [], 4 ;...
        handles.FST_err_edit,        'off', [], 'off', [], 'off', [], 4 ;...
        handles.FST_err_txt,         'off', [], 'off', [], 'off', [], 4 ;...
        handles.showFDsp_checkbox,   'off', 0 , 'off', 0 , 'off', [], 5 ;...
        handles.Waveform_popup,      'off', 1 , 'off', 1 , 'off', [], 6 ;...
        handles.Waveform_txt,        'off', [], 'off', [], 'off', [], 6 ;...
        handles.Freq_edit,           'off', [], 'on' , [], 'off', [], 7 ;...
        handles.Freq_txt,            'off', [], 'on' , [], 'off', [], 7 ;...
        handles.NHT_radiobutton,     'on' , 0 , 'on' , 0 , 'off', [], 8 ;...
        handles.NHT_tol_edit,        'off', [], 'off', [], 'off', [], 8 ;...
        handles.NHT_tol_txt,         'off', [], 'off', [], 'off', [], 8 ;...
        handles.FHT_radiobutton,     'on' , 1 , 'on' , 1 , 'off', [], 8 ;...
        handles.FHT_err_edit,        'on' , [], 'on' , [], 'off', [], 8 ;...
        handles.FHT_err_txt,         'on' , [], 'on' , [], 'off', [], 8 ;...
        handles.FDspline_checkbox,   'off', 0 , 'off', 0 , 'off', [], 9 ;...
        handles.FDsp_NDEC_edit,      'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_NDEC_txt,       'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_Bmin_Bmax_txt,  'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_Bmin_edit,      'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_Bmax_edit,      'off', [], 'off', [], 'off', [], 9 ;...
        handles.Quad_tol_edit,       'off', [], 'on' , [], 'off', [], 10;...
        handles.Quad_tol_txt,        'off', [], 'on' , [], 'off', [], 10;...
        handles.Rspline_checkbox,    'off', 0 , 'on' , 1 , 'off', [], 11;...
        handles.Rsp_NDEC_txt,        'off', [], 'on' , [], 'off', [], 11;...
        handles.Rsp_NDEC_edit,       'off', [], 'on' , [], 'off', [], 11;...
        handles.Seg_tol_edit,        'off', [], 'off', [], 'off', [], 12;... 
        handles.Seg_tol_txt,         'off', [], 'off', [], 'off', [], 12;... 
        handles.Max_seg_edit,        'off', [], 'off', [], 'off', [], 12;...
        handles.Max_seg_txt,         'off', [], 'off', [], 'off', [], 12};%    
Settings(2).includeInComparison = find(~ismember(                       ...
    [1:size(Settings(2).Defaults,1)], [9 10 12 13 20 21 23 24 36:39]));

Settings(3).Config  = {'Wenner';
                       'Schlumberger'};
Settings(3).defaultDomain = 1; % = DC
Settings(3).Defaults = {...        %  DC	 FD	    TD	      group
        handles.DC_radiobutton,      'on' , 1 , 'off', [], 'off', [], 1 ;...
        handles.Freqdom_radiobutton, 'off', 0 , 'off', [], 'off', [], 1 ;...
        handles.Timedom_radiobutton, 'off', 0 , 'off', [], 'off', [], 1 ;...
        handles.Full_sol_radiobutton,'off', 0 , 'off', [], 'off', [], 2 ;...
        handles.Quasi_radiobutton,   'off', 0 , 'off', [], 'off', [], 2 ;...
        handles.Times_edit,          'off', [], 'off', [], 'off', [], 3 ;...
        handles.Times_txt,           'off', [], 'off', [], 'off', [], 3 ;...
        handles.NST_radiobutton,     'off', 0 , 'off', [], 'off', [], 4 ;...
        handles.NST_tol_edit,        'off', [], 'off', [], 'off', [], 4 ;...
        handles.NST_tol_txt,         'off', [], 'off', [], 'off', [], 4 ;...
        handles.FST_radiobutton      'off', 0 , 'off', [], 'off', [], 4 ;...
        handles.FST_err_edit,        'off', [], 'off', [], 'off', [], 4 ;...
        handles.FST_err_txt,         'off', [], 'off', [], 'off', [], 4 ;...
        handles.showFDsp_checkbox,   'off', 0 , 'off', [], 'off', [], 5 ;...
        handles.Waveform_popup,      'off', 1 , 'off', [], 'off', [], 6 ;...
        handles.Waveform_txt,        'off', [], 'off', [], 'off', [], 6 ;...
        handles.Freq_edit,           'off', [], 'off', [], 'off', [], 7 ;...
        handles.Freq_txt,            'off', [], 'off', [], 'off', [], 7 ;...
        handles.NHT_radiobutton,     'on' , 0 , 'off', [], 'off', [], 8 ;...
        handles.NHT_tol_edit,        'off', [], 'off', [], 'off', [], 8 ;...
        handles.NHT_tol_txt,         'off', [], 'off', [], 'off', [], 8 ;...
        handles.FHT_radiobutton,     'on' , 1 , 'off', [], 'off', [], 8 ;...
        handles.FHT_err_edit,        'on' , [], 'off', [], 'off', [], 8 ;...
        handles.FHT_err_txt,         'on' , [], 'off', [], 'off', [], 8 ;...
        handles.FDspline_checkbox,   'off', 0 , 'off', [], 'off', [], 9 ;...
        handles.FDsp_NDEC_edit,      'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_NDEC_txt,       'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_Bmin_Bmax_txt,  'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_Bmin_edit,      'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_Bmax_edit,      'off', [], 'off', [], 'off', [], 9 ;...
        handles.Quad_tol_edit,       'off', [], 'off', [], 'off', [], 10;...
        handles.Quad_tol_txt,        'off', [], 'off', [], 'off', [], 10;...
        handles.Rspline_checkbox,    'off', 0 , 'off', [], 'off', [], 11;...
        handles.Rsp_NDEC_txt,        'off', [], 'off', [], 'off', [], 11;...
        handles.Rsp_NDEC_edit,       'off', [], 'off', [], 'off', [], 11;...
        handles.Seg_tol_edit,        'off', [], 'off', [], 'off', [], 12;... 
        handles.Seg_tol_txt,         'off', [], 'off', [], 'off', [], 12;... 
        handles.Max_seg_edit,        'off', [], 'off', [], 'off', [], 12;...
        handles.Max_seg_txt,         'off', [], 'off', [], 'off', [], 12};%    
Settings(3).includeInComparison = find(~ismember(                       ...
    [1:size(Settings(3).Defaults,1)], [9 10 12 13 20 21 23 24  36:39]));

Settings(4).Config  = {'HCP FDEM (HLEM)'};
Settings(4).defaultDomain = 2; % = Frequency Domain
Settings(4).Defaults = {...        %  DC	 FD	    TD	      group
        handles.DC_radiobutton,      'off', [], 'off', 0 , 'off', [], 1 ;...
        handles.Freqdom_radiobutton, 'off', [], 'on' , 1 , 'off', [], 1 ;...
        handles.Timedom_radiobutton, 'off', [], 'off', 0 , 'off', [], 1 ;...
        handles.Full_sol_radiobutton,'off', [], 'on' , 0 , 'off', [], 2 ;...
        handles.Quasi_radiobutton,   'off', [], 'on' , 1 , 'off', [], 2 ;...
        handles.Times_edit,          'off', [], 'off', [], 'off', [], 3 ;...
        handles.Times_txt,           'off', [], 'off', [], 'off', [], 3 ;...
        handles.NST_radiobutton,     'off', [], 'off', 0 , 'off', [], 4 ;...
        handles.NST_tol_edit,        'off', [], 'off', [], 'off', [], 4 ;...
        handles.NST_tol_txt,         'off', [], 'off', [], 'off', [], 4 ;...
        handles.FST_radiobutton      'off', [], 'off', 0 , 'off', [], 4 ;...
        handles.FST_err_edit,        'off', [], 'off', [], 'off', [], 4 ;...
        handles.FST_err_txt,         'off', [], 'off', [], 'off', [], 4 ;...
        handles.showFDsp_checkbox,   'off', [], 'off', 0 , 'off', [], 5 ;...
        handles.Waveform_popup,      'off', [], 'off', 1 , 'off', [], 6 ;...
        handles.Waveform_txt,        'off', [], 'off', [], 'off', [], 6 ;...
        handles.Freq_edit,           'off', [], 'on' , [], 'off', [], 7 ;...
        handles.Freq_txt,            'off', [], 'on' , [], 'off', [], 7 ;...
        handles.NHT_radiobutton,     'off', [], 'on' , 0 , 'off', [], 8 ;...
        handles.NHT_tol_edit,        'off', [], 'off', [], 'off', [], 8 ;...
        handles.NHT_tol_txt,         'off', [], 'off', [], 'off', [], 8 ;...
        handles.FHT_radiobutton,     'off', [], 'on' , 1 , 'off', [], 8 ;...
        handles.FHT_err_edit,        'off', [], 'on' , [], 'off', [], 8 ;...
        handles.FHT_err_txt,         'off', [], 'on' , [], 'off', [], 8 ;...
        handles.FDspline_checkbox,   'off', [], 'off', 0 , 'off', [], 9 ;...
        handles.FDsp_NDEC_edit,      'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_NDEC_txt,       'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_Bmin_Bmax_txt,  'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_Bmin_edit,      'off', [], 'off', [], 'off', [], 9 ;...
        handles.FDsp_Bmax_edit,      'off', [], 'off', [], 'off', [], 9 ;...
        handles.Quad_tol_edit,       'off', [], 'off', [], 'off', [], 10;...
        handles.Quad_tol_txt,        'off', [], 'off', [], 'off', [], 10;...
        handles.Rspline_checkbox,    'off', [], 'off', 0 , 'off', [], 11;...
        handles.Rsp_NDEC_txt,        'off', [], 'off', [], 'off', [], 11;...
        handles.Rsp_NDEC_edit,       'off', [], 'off', [], 'off', [], 11;...
        handles.Seg_tol_edit,        'off', [], 'off', [], 'off', [], 12;... 
        handles.Seg_tol_txt,         'off', [], 'off', [], 'off', [], 12;... 
        handles.Max_seg_edit,        'off', [], 'off', [], 'off', [], 12;...
        handles.Max_seg_txt,         'off', [], 'off', [], 'off', [], 12};%    
Settings(4).includeInComparison = find(~ismember(                       ...
    [1:size(Settings(4).Defaults,1)], [9 10 12 13 20 21 23 24  36:39]));

handles.Settings = Settings;

% Define context menu of the Calculate-button
cmenu = uicontextmenu('parent', handles.figure1);
set(handles.Calculate_button, 'UIContextMenu', cmenu);
cb = ['compute(''calculate_batch'',gcbo,[],guidata(gcbo));'];
handles.calc_batch_menu = uimenu(cmenu, 'Label', ['Batch calculation'] , 'Callback', cb, 'visible','off');

% Define context menu of the batch-button
cmenu = uicontextmenu('parent', handles.figure1);
set(handles.Batch_button, 'UIContextMenu', cmenu);
cb = ['compute(''load_batch'',gcbo,[],guidata(gcbo));'];
handles.load_batch_menu = uimenu(cmenu, 'Label', ['Load batch list'] , 'Callback', cb, 'visible','on');
cb = ['compute(''save_batch'',gcbo,[],guidata(gcbo));'];
handles.save_batch_menu = uimenu(cmenu, 'Label', ['Save batch list'] , 'Callback', cb, 'visible','off');
cb = ['compute(''clear_batch'',gcbo,[],guidata(gcbo));'];
handles.clear_batch_menu = uimenu(cmenu, 'Label', ['Clear batch list'] , 'Callback', cb, 'visible','off');

modelwin_han = guidata(handles.modelwin);
menu1 = [handles.calc_batch_menu; handles.save_batch_menu; handles.clear_batch_menu];
if isfield(modelwin_han, 'batch') && ~isempty(modelwin_han.batch)
    set(menu1, 'visible', 'on');
else
    set(menu1, 'visible', 'off');
end

guidata(hObject, handles);


% -------------------------------------------------------------------------
function setup_gui(hObject, eventdata, handles, domain)

modelwin_han = guidata(handles.modelwin);

k = 1;
while ~any(ismember(handles.Settings(k).Config,                         ...
        modelwin_han.config.type)) &&                                   ...
        ~(k>length(handles.Settings))
    k = k+1;
end % while, k now should hold index into settings for the current config

if k>length(handles.Settings)
    error('test','Invalid configuration!')
end % if

PropName = {'Enable', 'Value'};
id = [];

switch upper(domain)
    case {'DEFAULT'}
        % Choose either default domain, or current domain if one is
        % chosen and if it is allowed for the new config.
        if strcmpi(get(handles.figure1, 'visible'), 'on')
            currentDomain = cell(3,2);
            currentDomain = get(cell2mat(                               ...
                handles.Settings(k).Defaults(1:3,1)), PropName);
            m = 1;
            while currentDomain{m,2}~=1
                m = m+1;
            end
            if strcmpi(handles.Settings(k).Defaults(m,                  ...
                    handles.Settings(k).defaultDomain*2),               ...
                    currentDomain(m,1))
                id = [m*2,m*2+1];
            end % if
        else
            % If this is first run, make sure all UIcontrols are turned
            % off to ensure correct setting of default values below
            set([handles.Settings(1).Defaults{:,1}], 'enable', 'off');
        end % if
        if isempty(id)    
            id = [handles.Settings(k).defaultDomain.*2                  ...
                    handles.Settings(k).defaultDomain.*2+1];
        end % if
    case {'DC'}
        id = [2:3];
    case {'FD'}
        id = [4:5];
    case {'TD'}
        id = [6:7];
    otherwise
        error('Compute:setup_gui:domainerr','Invalid calculation domain!');
end % switch

% The following code makes sure we only change the values of the
% UIcontrols that must actually change 'enable' property. Thus even 
% though the default value of NHT_radiobutton is 0, if you have 
% previously set it to 1, and the 'enable' property does not change, 
% it should still be 1. 
CurrentState = cell(size(handles.Settings(k).Defaults,1),2);
CurrentState = get(cell2mat(handles.Settings(k).Defaults(:,1)), PropName);
GroupNos = cell2mat(handles.Settings(k).Defaults(                       ...
    handles.Settings(k).includeInComparison,8));
% find UIcontrols with a change in the 'Enable' state
mustChangeStateGroupNo = unique(GroupNos(                               ...
    ~strcmpi(handles.Settings(k).Defaults(                              ...
    handles.Settings(k).includeInComparison,id(1)),                     ...
    CurrentState(handles.Settings(k).includeInComparison,1))));

if ~isempty(mustChangeStateGroupNo)
    % find index into handles.Settings.Defaults
    GroupNos = cell2mat(handles.Settings(k).Defaults(:,8));
    mustChangeStateId = find(ismember(GroupNos, mustChangeStateGroupNo));
    % Set UICONTROLS settings
    set(cell2mat(handles.Settings(k).Defaults(mustChangeStateId,1)),    ...
        PropName(1:2), handles.Settings(k).Defaults(mustChangeStateId,id));
end

storeAllSettings(hObject, [], handles);


% -------------------------------------------------------------------------
function storeAllSettings(hObject, eventdata, handles)

% This function stores all relevant (selected) values and settings of the
% gui in the cparams structure, ready to pass on to one of the modelling
% routines.

if get(handles.DC_radiobutton, 'value') == 1
    cparams.domain = 'DC';
elseif get(handles.Freqdom_radiobutton, 'value') == 1
    cparams.domain = 'FD';
else
    cparams.domain = 'TD';
end % if

if get(handles.Full_sol_radiobutton, 'value') == 1
    cparams.calc_type = 'Full';
else
    cparams.calc_type = 'Quasi';
end % if

if cparams.domain == 'TD'
    cparams.times = get(handles.Times_edit, 'String');
    [tmp, isOk] = str2num(cparams.times);
    if ~isOk
        cparams.times = handles.DefaultStrings{                         ...
                [handles.DefaultStrings{:,1}]==handles.Times_edit,2};
    end % if
    if get(handles.NST_radiobutton, 'value') == 1
        cparams.FtoTtype = 'NST';
        [cparams.NST_tol, isOk] = str2num(get(handles.NST_tol_edit,     ...
            'String'));
        if ~isOk
            cparams.NST_tol = handles.DefaultStrings{                   ...
                    [handles.DefaultStrings{:,1}]==handles.NST_tol_edit,2};
        end % if
    else
        cparams.FtoTtype = 'FST';
        [cparams.FST_err, isOk] = str2num(get(handles.FST_err_edit,     ...
            'String'));
        if ~isOk
            cparams.FST_err = handles.DefaultStrings{                   ...
                    [handles.DefaultStrings{:,1}]==handles.FST_err_edit,2};
        end % if
    end % if
    
    cparams.showFDsp = get(handles.showFDsp_checkbox, 'Value');
    
    waveform = get(handles.Waveform_popup, 'String');
    cparams.waveform = waveform{get(handles.Waveform_popup, 'Value')};

elseif cparams.domain == 'FD'
    % Frequency domain specific parameters
    cparams.freq = get(handles.Freq_edit, 'String');
    [tmp, isOk] = str2num(cparams.freq);
    if ~isOk
        cparams.freq = handles.DefaultStrings{                         ...
                [handles.DefaultStrings{:,1}]==handles.Freq_edit,2};
    end % if
end

% Hankel transform parameters
if get(handles.NHT_radiobutton, 'value') == 1
    cparams.hank_type = 'NHT';
    [cparams.NHT_tol, isOk] = str2num(get(handles.NHT_tol_edit,         ...
        'String'));
    if ~isOk
        cparams.NHT_tol = handles.DefaultStrings{                       ...
                [handles.DefaultStrings{:,1}]==handles.NHT_tol_edit,2};
    end % if
else
    cparams.hank_type = 'FHT';
    [cparams.FHT_err, isOk] = str2num(get(handles.FHT_err_edit,         ...
        'String'));
    if ~isOk
        cparams.FHT_err = handles.DefaultStrings{                       ...
                [handles.DefaultStrings{:,1}]==handles.FHT_err_edit,2};
    end % if
end % if

cparams.FDspline = get(handles.FDspline_checkbox, 'Value');
if cparams.FDspline
    cparams.FDsp_NDEC = str2double(get(handles.FDsp_NDEC_edit, 'String'));
    cparams.FDsp_Bmin = str2double(get(handles.FDsp_Bmin_edit, 'String'));
    cparams.FDsp_Bmax = str2double(get(handles.FDsp_Bmax_edit, 'String'));
end % if

if strcmpi(get(handles.Quad_tol_edit, 'Enable'), 'On')
    cparams.Quad_tol = str2double(get(handles.Quad_tol_edit, 'String'));
    cparams.Rspline = get(handles.Rspline_checkbox, 'Value');
    if cparams.Rspline
        cparams.Rsp_NDEC = str2double(get(                              ...
            handles.Rsp_NDEC_edit, 'String'));
    end % if
end % if

if isfield(cparams, 'hank_type') && strcmpi(cparams.hank_type, 'NHT') ||...
        isfield(cparams, 'FtoTtype') && strcmpi(cparams.FtoTtype, 'NST')
    cparams.Seg_tol = str2double(get(handles.Seg_tol_edit, 'String'));
    cparams.Max_seg = str2double(get(handles.Max_seg_edit, 'String'));
end

modelwin_han = guidata(handles.modelwin);
modelwin_han.cparams = [];
modelwin_han.cparams = cparams;
guidata(handles.modelwin, modelwin_han);

% -------------------------------------------------------------------------
function loadAllSettings(hObject, eventdata, handles,cparams)

% This function loads all relevant (selected) values and settings from the cparams structure
% to the gui.

domain_set = [handles.DC_radiobutton,handles.Freqdom_radiobutton,handles.Timedom_radiobutton];
switch cparams.domain
    case {'DC'}
        set(domain_set, {'value'}, {1;0;0});
    case {'FD'}
        set(domain_set, {'value'}, {0;1;0});
    case {'TD'}
        set(domain_set, {'value'}, {0;0;1});
end

if isfield(cparams,'calc_type')
    calc_set = [handles.Quasi_radiobutton,handles.Full_sol_radiobutton];
    switch cparams.calc_type
        case {'Quasi'}
            set(calc_set, {'value'}, {1;0});
        case {'Full'}
            set(calc_set, {'value'}, {0;1});
    end
end

if isfield(cparams,'times')
    set(handles.Times_edit, 'String', mat2str(cparams.times));
end

if isfield(cparams,'FtoTtype')
    FtoT_set = [handles.NST_radiobutton,handles.FST_radiobutton];
    switch cparams.FtoTtype
        case {'NST'}
            set(FtoT_set, {'value'}, {1;0});
            set(handles.NST_tol_edit, 'string', sprintf('%e',cparams.NST_tol));
        case {'FST'}
            set(FtoT_set, {'value'}, {0;1});
            set(handles.FST_err_edit, 'string', sprintf('%e',cparams.FST_err));
    end
end

if isfield(cparams,'showFDsp')
    set(handles.showFDsp_checkbox, 'Value', cparams.showFDsp);
end
    
if isfield(cparams,'waveform')
    waveform = get(handles.Waveform_popup, 'String');
    set(handles.Waveform_popup, 'Value', find(ismember(cparams.waveform, waveform)));
end
    
if isfield(cparams,'freq')
    set(handles.Freq_edit, 'String', mat2str(cparams.freq));
end

if isfield(cparams,'hank_type')
    hank_set = [handles.NHT_radiobutton,handles.FHT_radiobutton];
    switch cparams.hank_type
        case {'NHT'}
            set(hank_set, {'value'}, {1;0});
            set(handles.NHT_tol_edit, 'string', sprintf('%e',cparams.NHT_tol));
        case {'FST'}
            set(hank_set, {'value'}, {0;1});
            set(handles.FHT_err_edit, 'string', sprintf('%e',cparams.FHT_err));
    end
end

if isfield(cparams,'FDspline')
    set(handles.FDspline_checkbox, 'Value', cparams.FDspline);
    if cparams.FDspline
        set(handles.FDsp_NDEC_edit, 'String', sprintf('%d', cparams.FDsp_NDEC));
        set(handles.FDsp_Bmin_edit, 'String', sprintf('%d', cparams.FDsp_Bmin));
        set(handles.FDsp_Bmax_edit, 'String', sprintf('%d', cparams.FDsp_Bmax));
    end % if
end

if isfield(cparams,'Quad_tol')
    set(handles.Quad_tol_edit, 'String', sprintf('%e',cparams.Quad_tol));
end
if isfield(cparams,'Rspline')
    set(handles.Rspline_checkbox, 'Value', cparams.Rspline);
    if cparams.Rspline
        set(handles.Rsp_NDEC_edit, 'String', sprintf('%d', cparams.Rsp_NDEC));
    end % if
end

if isfield(cparams,'Seg_tol')
    set(handles.Seg_tol_edit, 'String', sprintf('%e', cparams.Seg_tol));
end
if isfield(cparams,'Max_seg')
    set(handles.Max_seg_edit, 'String', sprintf('%d', cparams.Max_seg));
end

storeAllSettings(hObject, [], handles);


% *************************************************************************
% *                                                                       *
% *   Callbacks for Calculation UICONTROLS                                *
% *                                                                       *
% *************************************************************************

% -------------------------------------------------------------------------
function DC_radiobutton_Callback(hObject, eventdata, handles)

set(handles.DC_radiobutton, 'Value', 1);
set(handles.Freqdom_radiobutton, 'Value', 0);
set(handles.Timedom_radiobutton, 'Value', 0);
setup_gui(hObject, eventdata, handles, 'DC');

% -------------------------------------------------------------------------
function Freqdom_radiobutton_Callback(hObject, eventdata, handles)

set(handles.DC_radiobutton, 'Value', 0);
set(handles.Freqdom_radiobutton, 'Value', 1);
set(handles.Timedom_radiobutton, 'Value', 0);
setup_gui(hObject, eventdata, handles, 'FD');

% -------------------------------------------------------------------------
function Timedom_radiobutton_Callback(hObject, eventdata, handles)

set(handles.DC_radiobutton, 'Value', 0);
set(handles.Freqdom_radiobutton, 'Value', 0);
set(handles.Timedom_radiobutton, 'Value', 1);
setup_gui(hObject, eventdata, handles, 'TD');

% -------------------------------------------------------------------------
function Full_sol_radiobutton_Callback(hObject, eventdata, handles)

set(handles.Quasi_radiobutton, 'Value', 0);
set(handles.Full_sol_radiobutton, 'Value', 1);

% -------------------------------------------------------------------------
function Quasi_radiobutton_Callback(hObject, eventdata, handles)

set(handles.Quasi_radiobutton, 'Value', 1);
set(handles.Full_sol_radiobutton, 'Value', 0);


% *************************************************************************
% *                                                                       *
% *   Callbacks for Frequency parameters UICONTROLS                       *
% *                                                                       *
% *************************************************************************

% -------------------------------------------------------------------------
function Freq_edit_Callback(hObject, eventdata, handles)

[num, isOk] = str2num(get(handles.Freq_edit, 'string'));

if ~isOk
    modelwin_han = guidata(handles.modelwin);
    set(handles.Freq_edit, 'string', modelwin_han.cparams.freq);
end


% -------------------------------------------------------------------------
function NHT_radiobutton_Callback(hObject, eventdata, handles)

set(handles.NHT_radiobutton, 'Value', 1);
set(handles.FHT_radiobutton, 'Value', 0);
set(handles.NHT_tol_edit,    'Enable', 'on' );
set(handles.NHT_tol_txt,     'Enable', 'on' );
set(handles.Max_seg_edit,    'Enable', 'on' );
set(handles.Max_seg_txt,     'Enable', 'on' );
set(handles.Seg_tol_edit,    'Enable', 'on' );
set(handles.Seg_tol_txt,     'Enable', 'on' );
set(handles.FHT_err_edit,    'Enable', 'off');
set(handles.FHT_err_txt,     'Enable', 'off');


% -------------------------------------------------------------------------
function FHT_radiobutton_Callback(hObject, eventdata, handles)

set(handles.NHT_radiobutton, 'Value', 0);
set(handles.FHT_radiobutton, 'Value', 1);
set(handles.NHT_tol_edit,    'Enable', 'off');
set(handles.NHT_tol_txt,     'Enable', 'off');
if get(handles.FST_radiobutton, 'Value') ||                             ...
        strcmpi(get(handles.FST_radiobutton, 'Enable'), 'off')
    set(handles.Max_seg_edit,    'Enable', 'off');
    set(handles.Max_seg_txt,     'Enable', 'off');
    set(handles.Seg_tol_edit,    'Enable', 'off');
    set(handles.Seg_tol_txt,     'Enable', 'off');
end
set(handles.FHT_err_edit,    'Enable', 'on' );
set(handles.FHT_err_txt,     'Enable', 'on' );


% -------------------------------------------------------------------------
function FDspline_checkbox_Callback(hObject, eventdata, handles)

if get(handles.FDspline_checkbox, 'Value')
    set(handles.showFDsp_checkbox, 'Enable', 'on' , 'Value', 0);
    set(handles.FDsp_NDEC_edit,        'Enable', 'on' );
    set(handles.FDsp_NDEC_txt,         'Enable', 'on' );
    set(handles.FDsp_Bmin_Bmax_txt,    'Enable', 'on' );
    set(handles.FDsp_Bmin_edit,        'Enable', 'on' );
    set(handles.FDsp_Bmax_edit,        'Enable', 'on' );
else
    set(handles.showFDsp_checkbox, 'Enable', 'off', 'Value', 0);
    set(handles.FDsp_NDEC_edit,        'Enable', 'off');
    set(handles.FDsp_NDEC_txt,         'Enable', 'off');
    set(handles.FDsp_Bmin_Bmax_txt,    'Enable', 'off');
    set(handles.FDsp_Bmin_edit,        'Enable', 'off');
    set(handles.FDsp_Bmax_edit,        'Enable', 'off');
end


% *************************************************************************
% *                                                                       *
% *   Callbacks for Time parameters UICONTROLS                            *
% *                                                                       *
% *************************************************************************

% -------------------------------------------------------------------------
function Times_edit_Callback(hObject, eventdata, handles)

[num, isOk] = str2num(get(handles.Times_edit, 'string'));

if ~isOk
    modelwin_han = guidata(handles.modelwin);
    set(handles.Times_edit, 'string', modelwin_han.cparams.times);
end


% -------------------------------------------------------------------------
function numberinput_Callback(hObject, eventdata, handles)

num = str2double(get(hObject, 'string'));

if isnan(num)
    set(hObject, 'string', handles.DefaultStrings{                      ...
            [handles.DefaultStrings{:,1}]==hObject,2});
end


% -------------------------------------------------------------------------
function NST_radiobutton_Callback(hObject, eventdata, handles)

set(handles.NST_radiobutton, 'Value', 1);
set(handles.FST_radiobutton, 'Value', 0);
set(handles.FST_err_edit,    'Enable', 'off');
set(handles.FST_err_txt,     'Enable', 'off');
set(handles.NST_tol_edit,    'Enable', 'on' );
set(handles.NST_tol_txt,     'Enable', 'on' );
set(handles.Max_seg_edit,    'Enable', 'on' );
set(handles.Max_seg_txt,     'Enable', 'on' );
set(handles.Seg_tol_edit,    'Enable', 'on' );
set(handles.Seg_tol_txt,     'Enable', 'on' );

modelwin_han = guidata(handles.modelwin);
modelwin_han.cparams.FtoTtype = 'NST';
guidata(handles.modelwin,modelwin_han);


% -------------------------------------------------------------------------
function FST_radiobutton_Callback(hObject, eventdata, handles)

set(handles.NST_radiobutton, 'Value', 0);
set(handles.FST_radiobutton, 'Value', 1);
set(handles.NST_tol_edit,    'Enable', 'off');
set(handles.NST_tol_txt,     'Enable', 'off');
if get(handles.FHT_radiobutton, 'Value') ||                             ...
        strcmpi(get(handles.FHT_radiobutton, 'Enable'), 'off')
    set(handles.Max_seg_edit,    'Enable', 'off');
    set(handles.Max_seg_txt,     'Enable', 'off');
    set(handles.Seg_tol_edit,    'Enable', 'off');
    set(handles.Seg_tol_txt,     'Enable', 'off');
end
set(handles.FST_err_edit,    'Enable', 'on' );
set(handles.FST_err_txt,     'Enable', 'on' );


% -------------------------------------------------------------------------
function showFDsp_checkbox_Callback(hObject, eventdata, handles)


% -------------------------------------------------------------------------
function Waveform_popup_Callback(hObject, eventdata, handles)
set(hObject, 'Value', 1);
disp('Waveform popup not implemented yet...');


% *************************************************************************
% *                                                                       *
% *   Callbacks for Dipole integrations UICONTROLS                        *
% *                                                                       *
% *************************************************************************

% -------------------------------------------------------------------------
function Rspline_checkbox_Callback(hObject, eventdata, handles)

if get(handles.Rspline_checkbox, 'Value')
    set(handles.Rsp_NDEC_edit,   'Enable', 'on' );
    set(handles.Rsp_NDEC_txt,    'Enable', 'on' );
else
    set(handles.Rsp_NDEC_edit,   'Enable', 'off');
    set(handles.Rsp_NDEC_txt,    'Enable', 'off');
end


% *************************************************************************
% *                                                                       *
% *   Callbacks for Pushbutton UICONTROLS                                 *
% *                                                                       *
% *************************************************************************

% --------------------------------------------------------------------
function load_batch(hObject, eventdata, handles)

modelwin_han = guidata(handles.modelwin);
if isfield(modelwin_han,'batch') && ~isempty(modelwin_han.batch)
    button = questdlg('Current list will be overwritten!',...
        'Clear list?','OK','Cancel','Cancel');
    if strcmp(button,'Cancel')
        return;
    end
end

[filename, pathname] = uigetfile(                                       ...
    {'*.mat','MAT-files (*.mat)';                                       ...
        '*.*',  'All Files (*.*)'},                                     ...
    'Save batch list as...', '*.mat');
if isequal(filename,0), return; end
file = fullfile(pathname, filename);
if exist(file,'file')
    batchlist = load(file,'batchlist');
    if ~isfield(batchlist,'batchlist'), return; end
    batchlist = batchlist.batchlist;
    names = {'cparams';'config';'layers'};
    if ~all(ismember(names, fieldnames(batchlist)))
        disp('Bad file!');
        disp('Aborting...');
        return;
    else
        modelwin_han.batch = batchlist;
    end
end

menu1 = [handles.calc_batch_menu; handles.save_batch_menu; handles.clear_batch_menu];
if ~isempty(modelwin_han.batch)
    set(menu1, 'visible', 'on');
else
    set(menu1, 'visible', 'off');
end
guidata(handles.modelwin, modelwin_han);


% --------------------------------------------------------------------
function save_batch(hObject, eventdata, handles)

[filename, pathname] = uiputfile(                                       ...
            {'*.mat','MAT-files (*.mat)';                                       ...
                '*.*',  'All Files (*.*)'},                                     ...
            'Save batch list as...', 'batch_list.mat');

if isequal(filename,0), return; end
[tmp1,filename,ext] = fileparts(filename);
if isempty(ext), ext = '.mat'; end
modelwin_han = guidata(handles.modelwin);
batchlist = modelwin_han.batch; 
save(fullfile(pathname, [filename ext] ),'batchlist','-mat');     

% --------------------------------------------------------------------
function clear_batch(hObject, eventdata, handles)

button = questdlg('Clear batch list?',...
    'Clear list?','Yes','No','No');
if strcmp(button,'Yes')
    modelwin_han = guidata(handles.modelwin);
    modelwin_han.batch = [];
    guidata(handles.modelwin, modelwin_han);
    
    menu1 = [handles.calc_batch_menu; handles.save_batch_menu; handles.clear_batch_menu];
    set(menu1, 'visible', 'off');
    
elseif strcmp(button,'No')
    return;
end    


% --------------------------------------------------------------------
function Batch_button_Callback(hObject, eventdata, handles)

storeAllSettings(hObject, [], handles);
modelwin_han = guidata(handles.modelwin);
model.cparams = modelwin_han.cparams;
model.config = rmfield(modelwin_han.config,{'plot_handle'});
model.layers = rmfield(modelwin_han.layers,{'label_handle'});

if ~isfield(modelwin_han, 'batch')
    modelwin_han.batch = [];
end

[model, batmodels] = batch(model, modelwin_han.batch);

if ~isempty(model) && ~isempty(batmodels)
    modelwin_han.batch = batmodels;
    if isfield(model, 'name')
        modelwin_han.name = model.name;
    end
    
    menu1 = [handles.calc_batch_menu; handles.save_batch_menu; handles.clear_batch_menu];
    if ~isempty(modelwin_han.batch)
        set(menu1, 'visible', 'on');
    else
        set(menu1, 'visible', 'off');
    end
    
    guidata(handles.modelwin, modelwin_han);
    guidata(hObject, handles);
end
 

% --------------------------------------------------------------------
function pushbutton13_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Time_contextmenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Waveform_cmenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Load_cmenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Gates_cmenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function Save_cmenu_Callback(hObject, eventdata, handles)


% -------------------------------------------------------------------------
function Freq_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function NHT_tol_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function FHT_err_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function FDsp_Bmin_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function FDsp_Bmax_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function FDsp_NDEC_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function Times_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function NST_tol_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function Quad_tol_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function Rsp_NDEC_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function Waveform_popup_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% -------------------------------------------------------------------------
function FST_err_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function Seg_tol_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes during object creation, after setting all properties.
function Max_seg_edit_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% *************************************************************************
% *                                                                       *
% *   Functions for storing and plotting the results                      *
% *                                                                       *
% *************************************************************************

% function plot_tem(v, config, cparams, fig)
% 
% if getappdata(0, 'debug')
%     disp('comput2.m:Plot_TEM: Debug before plot...');
% end
% 
% figure(fig, 'visible', 'on');
% 
% if strcmp(cparams.domain,'TD')
%     neg_loglog(cparams.times,v,[],'Markersize',5);
%     ylabel('V/A');
%     % rhoa = (config.TxR.^(4/3).*config.RxA.^(2/3).*(4e-7*pi).^(5/3))./ ...
%     %     (20.^(2/3).*pi.^(1/3).*cparams.times.^(5/3).*v.^(2/3));
%     % loglog(cparams.times,rhoa,'-xb','Markersize',5);
%     % ylabel('\rho_a (\Omega m)');
%     xlabel('Time (sec)');
% else
%     subplot(2,1,1);
%     semilogx(cparams.freq, real(v), '-xb', 'Markersize',5);
%     ylabel('Real(H_z)');
%     subplot(2,1,2);
%     semilogx(cparams.freq, imag(v), '-xb', 'Markersize',5);
%     ylabel('Imag(H_z)');
%     xlabel('Frequency (Hz)');
% end


% --------------------------------------------------------------------
% function plot_hcp(result, config, cparams, fig)
% 
% if min(size(result)) == 1
%     figure(fig, 'visible', 'on');                         
%     for k = 1:length(result)
%         subplot(2,1,1);
%         semilogx(cparams.freq, real(result(k).H./result(k).H_prim));
%         hold on;
%         subplot(2,1,2);
%         semilogx(cparams.freq, imag(result(k).H./result(k).H_prim));
%         hold on;
%     end % k
%     subplot(2,1,1);
%     ylabel('Real(Z/Z_0)');
%     hold off;
%     subplot(2,1,2);
%     ylabel('Imag(Z/Z_0)');
%     xlabel('Freq. (Hz)');
%     hold off;         
% else
%     disp('I don'' know how to plot these results!');
% end % if


% --------------------------------------------------------------------
% function plot_emgsa(result, config, cparams, fig)
% 
% if min(size(result)) == 1
%     figure(fig);
%     set(fig, 'visible', 'on');                        
%     for k = 1:length(result)
%         if ~isfield(result(k), 'G_factor')        
%             result(k).G_factor = 0;
%             if isfield(result(k), 'Aspac')
%                 result(k).G_factor = pi.*result(k).Aspac                ...
%                     .*(result(k).Nspac)                                 ...
%                     .*(result(k).Nspac+1)                               ...
%                     .*(result(k).Nspac+2);
%             end
%             if isempty(result(k).G_factor) || result(k).G_factor == 0
%                 G_factor = 1;
%             end % if
%         end
%         subplot(2,1,1);
% %        halldat = - 6.7e-011*cparams.freq.^3 +                          ...
% %            1.1e-006*cparams.freq.^2 - 0.01*cparams.freq - 1;
%         semilogx(cparams.freq, abs(result(k).Z.*result(k).G_factor));
%         hold on;
%         subplot(2,1,2);
%         semilogx(cparams.freq,                                          ...
%             -angle(result(k).Z).*1000);%+halldat.');
%         hold on;
%     end % k
%     subplot(2,1,1);
%     ylabel('\rho_a (\Omega\cdot m)');
%     hold off;
%     subplot(2,1,2);
%     ylabel('neg. phase (mrad)');
%     xlabel('Freq. (Hz)');
%     hold off;         
% else
%     % if you want results in a 3D matrix for D-D data
%     % Z = reshape([result.Z],length(result(1).Z),  ...
%     %    size(result,1),size(result,2));
%     disp('I don'' know how to plot these results!');
% end % if


% --------------------------------------------------------------------
% function plot_dcgsa(result, config, cparams, fig)
% 
% if min(size(result)) == 1 && max(size(result)) > 1
%     figure(fig, 'visible', 'on');
%     yval = squeeze(reshape([result.Z], size(result,1), size(result,2)));
%     G_factor = squeeze(reshape([result.G_factor], size(result,1),       ...
%         size(result,2)));
%     
%     for k = 1:length(yval)
%         if G_factor(k) == 0
%             G_factor(k) = 1;
%         end % if
%         switch config.type
%             case {'Dipole-Dipole', '*Capacitance*'}
%                 if size(result,1) > 1
%                     xval = squeeze(reshape(         ...
%                         [result.Aspac],             ...
%                         size(result,1), 1));
%                     loglog(xval, yval.*G_factor);
%                 else
%                     xval = squeeze(reshape(         ...
%                         [result.Nspac],             ...
%                         1, size(result,2)));
%                     semilogy(xval, yval.*G_factor);
%                 end % if
%             case {'Wenner'}
%                 xval = squeeze(reshape(             ...
%                     [result.Aspac],                 ...
%                     size(result,1), 1));
%                 semilogy(xval, yval.*G_factor);
%             case {'Schlumberger'}
%                 if size(result,1) > 1
%                     xval = squeeze(reshape(         ...
%                         [result.OA],                ...
%                         size(result,1), 1));
%                     loglog(xval, yval.*G_factor);                                                                            
%                 else
%                     xval = squeeze(reshape(         ...
%                         [result.OM],                ...
%                         1, size(result,2)));
%                     semilogy(xval, yval.*G_factor);
%                 end % if
%         end
%         hold on;
%     end % k
%     hold off;
%     ylabel('Apparent Resistivity (\Omegam)');
%     switch config.type
%         case {'Dipole-Dipole', '*Capacitance*'}
%             if size(result,1) > 1
%                 xlabel('A-spacing (m)');
%             else
%                 xlabel('N-spacing (m)');
%             end % if
%         case {'Wenner'}
%             xlabel('A-spacing (m)');
%         case {'Schlumberger'}
%             if size(result,1) > 1
%                 xlabel('OA (m)');                                                                            
%             else
%                 xlabel('OM (m)');
%             end % if
%     end % switch   
% end % if


% --------------------------------------------------------------------
function saveforward(model, varargin)

if getappdata(0, 'debug')
    disp('compute.m:saveforwad:debug:  Entered saveforward routine');        
end

if nargin == 2
    [pathname,filename,ext] = fileparts(varargin{1});
    filterindex = [];
    switch lower(ext)
        case {'','.mat'}
            filterindex = 1;
        case {'.dat'}
            filterindex = 2;
        otherwise
            filterindex = 3;
    end
else    
    
    if getappdata(0, 'compiled')
        % !!! for standalone
        [filename, pathname] = uiputfile(                                       ...
            {'*.mat','MAT-files (*.mat)';                                       ...
                '*.dat','ASCII file (*.dat)';                                   ...
                '*.*',  'All Files (*.*)'},                                     ...
            'Save result as...', 'model_result.mat');
        if isequal(filename,0), return; end
        
        filterindex = [];
        [tmp1,filename,ext] = fileparts(filename);
        switch lower(ext)
            case {'','.mat'}
                filterindex = 1;
            case {'.dat'}
                filterindex = 2;
            otherwise
                filterindex = 3;
        end
        
    else
        % !!! for matlab mode
        [filename, pathname, filterindex] = uiputfile(                          ...
            {'*.mat','MAT-files (*.mat)';                                       ...
                '*.dat','ASCII file (*.dat)';                                   ...
                '*.*',  'All Files (*.*)'},                                     ...
            'Save result as...', 'model_result.mat');
        if isequal(filename,0), return; end
        [tmp1,filename,ext] = fileparts(filename);
    end
end

% model.config = rmfield(handles.config,{'plot_handle'});
% model.layers = rmfield(handles.layers,{'label_handle'});
% if isfield(model.layers,'interface_handle')
%     model.layers = rmfield(model.layers,{'interface_handle',            ...
%             'prop_lab_handle'});
% end

if filterindex == 1 || filterindex == 3
    if isempty(ext), ext = '.mat'; end
    save(fullfile(pathname, [filename ext] ),'model','-mat');
elseif filterindex==2
    fid = fopen(fullfile(pathname, [filename ext]),'w');
    fprintf(fid,'Configuration:  %s\n',model.config.type);
    switch model.config.type
        case {'Dipole-Dipole', '*Capacitance*'}
            fprintf(fid,'A-spacing:      ');
            fprintf(fid,'%f\t',model.config.Aspac);
            fprintf(fid,'\nN-spacing:      ');
            fprintf(fid,'%f\t',model.config.Nspac);
        case {'Schlumberger'}
            fprintf(fid,'OA:             %f\n',model.config.OA);
            fprintf(fid,'OM:             %f\n',model.config.OM);
        case {'Wenner'}
            fprintf(fid,'A-spacing:      %f\n',model.config.a);
        case {'HCP FDEM (HLEM)'}
            fprintf(fid,'R-spacing:      %f\n',model.config.Rspac);
        case {'TEM Central Loop'}
            fprintf(fid,'Tx-side:        %f\n',model.config.TxS);
            fprintf(fid,'Rx-radius:      %f\n',model.config.RxR);
    end
    
    fprintf(fid,'\nModel:\n');            
    fprintf(fid,'Rho DC:         ');
    fprintf(fid,'%f\t',[model.layers.rho]);
    fprintf(fid,'\nThickness:      ');
    fprintf(fid,'%f\t',[model.layers.thickness]);
    fprintf(fid,'\ndepth_to_top:   ');
    fprintf(fid,'%f\t',[model.layers.depth_to_top]);
    
    fprintf(fid,'\n\nCalculation: %s\n',model.calc_type);
    fprintf(fid,'Transform type:  %s\n',model.hank_type);
    if strcmp(model.hank_type,'cconvol')
        fprintf(fid,'cconvol tol:    %s\n',model.cconvolerr);
    else
        fprintf(fid,'Hank tol:       %s\n',model.Hank_tol);
    end
    
    fprintf(fid,'\nResult:\n');
    switch model.calc_type
        case{'DC'}
            switch model.config.type
                case {'Dipole-Dipole', '*Capacitance*'}                        
                    fprintf(fid,'Aspac\t\tNspac\t\tRhoA\n');                            
                    for k = 1:length(model.result)
                        fprintf(fid,'%f\t%f\t%f',               ...
                            model.result(k).Aspac,              ...
                            model.result(k).Nspac,              ...
                            model.result(k).Z);
                    end
                case {'Schlumberger'}
                    fprintf(fid,'OA\t\tOM\t\tRhoA\n');                            
                    for k = 1:length(model.result)
                        fprintf(fid,'%f\t%f\t%f',               ...
                            model.result(k).OA,                 ...
                            model.result(k).OM,                 ...
                            model.result(k).Z);
                    end
                case {'Wenner'}
                    for k = 1:length(model.result)
                        fprintf(fid,'%f\t%f',                   ...
                            model.result(k).Aspac,              ...
                            model.result(k).Z);
                    end
                case {'General Surface Array','*GSA capacitance*'}
                    for k = 1:length(model.result)
                        fprintf(fid,'%i\t%f',                   ...
                            k, model.result(k).Z);
                    end    
            end
        case {'FD'}
            switch model.config.type                
                case {'HCP FDEM (HLEM)'}
                    for k = 1:length(model.result)
                        fprintf(fid,'%f\t%f',                   ...
                            model.result(k).Rspac,              ...
                            model.result(k).Z);
                    end
            end
        otherwise
            disp('I don''t know how to save these data...!');
    end
    fclose(fid);
end



% --- Executes on button press in Reusefig_checkbox.
function Reusefig_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Reusefig_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Reusefig_checkbox


