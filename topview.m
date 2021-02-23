function varargout = topview(varargin)

% Usage:
%    config = topview(config);
%
% topview.m is part of the CR1Dmod forward modeling package. The function
% related to the topview.fig GUI, which allows the user to place electrodes
% and wires arbitrarily on the surface of the half-space.
% The user must ensure that no wires cross! There is no automatic 
% detection of this source of error implemented at this time.
%
% 'config' is an optional input structure with at least the fields:
%        type:    Electrode array ('Dipole-Dipole', 'Wenner',
%                                  'Schlumberger' or 'General Surface Array')
%        C1:      [x,y,z] coordinates of first current electrode
%        C2:      [x,y,z] coordinates of second current electrode
%        P1:      [x,y,z] coordinates of first potential electrode
%        P2:      [x,y,z] coordinates of second potential electrode
%
% topview.m returns a new config structure which in addition to the above 
% holds the fields Cwire and Pwire with coordinates of the additional points 
% defining the wires. The new structure also holds the fields Aspac and
% Nspace, which only has a physical significance in dipole-dipole and 
% wenner configrations.
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
                   'gui_OpeningFcn', @topview_OpeningFcn, ...
                   'gui_OutputFcn',  @topview_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --------------------------------------------------------------------
% --- Executes just before topview is made visible.
function topview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to topview (see VARARGIN)

% Choose default command line output for topview
handles.output = [];

% Update handles structure
guidata(hObject, handles);

set(handles.Add_button, 'enable', 'off');

handles.config.type = 'Dipole-Dipole';
handles.config.Aspac = 100;
handles.config.Nspac = 1;
handles.config.C1 = [-handles.config.Aspac-...
        0.5*handles.config.Aspac*handles.config.Nspac 0 0];  % transmitter dipole            C1  a  C2     P1  a  P2
handles.config.C2 = [handles.config.C1(1)+...
        handles.config.Aspac(1) 0 0];                        % transmitter dipole            o------o      o------o    
handles.config.P1 = [handles.config.C2(1)+...                 % receiver dipole                 Tx     n*a    Rx
        handles.config.Nspac(1)*handles.config.Aspac(1) 0 0];   
handles.config.P2 = [handles.config.P1(1)+...                 % receiver dipole
        handles.config.Aspac(1) 0 0];     
handles.config.Cwire = [];
handles.config.Pwire = [];

handles.topview.Cwire = line([handles.config.C1(1) handles.config.C2(1)],...
    [handles.config.C1(2) handles.config.C2(2)],'MarkerSize',6, ...
    'markeredgecolor','k','markerfacecolor','b','marker','none',...
    'linewidth',2,'color','k','ButtonDownFcn',{@apex_ButtonDownFcn, 'Cwire'});

handles.topview.Pwire = line([handles.config.P1(1) handles.config.P2(1)],...
    [handles.config.P1(2) handles.config.P2(2)],'MarkerSize',6, ...
    'markeredgecolor','k','markerfacecolor','b','marker','none',...
    'linewidth',2,'color','k','ButtonDownFcn',{@apex_ButtonDownFcn, 'Pwire'});

handles.topview.current_apex = {'Cwire', 0, []};

handles.topview.C1 = line(handles.config.C1(1), handles.config.C1(2),'MarkerSize',8, ...
    'markeredgecolor','k','markerfacecolor','r','marker','o',...
    'ButtonDownFcn','topview(''electrode_ButtonDownFcn'',gcbo,[],guidata(gcbo),''C1'');');
handles.topview.C2 = line(handles.config.C2(1), handles.config.C2(2),'MarkerSize',8, ...
    'markeredgecolor','k','markerfacecolor','r','marker','o',...
    'ButtonDownFcn','topview(''electrode_ButtonDownFcn'',gcbo,[],guidata(gcbo),''C2'');');
handles.topview.P1 = line(handles.config.P1(1), handles.config.P1(2),'MarkerSize',8, ...
    'markeredgecolor','k','markerfacecolor','b','marker','o',...
    'ButtonDownFcn','topview(''electrode_ButtonDownFcn'',gcbo,[],guidata(gcbo),''P1'');');
handles.topview.P2 = line(handles.config.P2(1), handles.config.P2(2),'MarkerSize',8, ...
    'markeredgecolor','k','markerfacecolor','b','marker','o',...
    'ButtonDownFcn','topview(''electrode_ButtonDownFcn'',gcbo,[],guidata(gcbo),''P2'');');

if length(varargin)==1
    names = fieldnames(varargin{1});
    for k = 1:length(names)
        handles.config.(names{k}) = varargin{1}.(names{k});
    end
    
    if ~isempty(handles.config.Cwire);
        pos = handles.config.Cwire;
        handles.config.Cwire = [];
        for k = 1:size(pos,1)
            handles = add_apex(handles,pos(k,1:2),'Cwire');
        end
    end
    
    if ~isempty(handles.config.Pwire);
        handles.topview.current_apex(1:2) = {'Pwire', 0};
        pos = handles.config.Pwire;
        handles.config.Pwire = [];
        for k = 1:size(pos,1)
            handles = add_apex(handles,pos(k,1:2),'Pwire');
        end
    end
end

set(handles.Config_popup,'string',handles.config.type);

handles = update_topview(handles);

set_axis_limmits(handles);

% Update handles structure
guidata(hObject, handles);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')


% UIWAIT makes topview wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = topview_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);



% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function C1_x_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to C1_x_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --------------------------------------------------------------------
function C1_x_edit_Callback(hObject, eventdata, handles)
% hObject    handle to C1_x_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of C1_x_edit as text
%        str2double(get(hObject,'String')) returns contents of C1_x_edit as a double
handles.config.C1(1) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function C1_y_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function C1_y_edit_Callback(hObject, eventdata, handles)
handles.config.C1(2) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function C1_z_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function C1_z_edit_Callback(hObject, eventdata, handles)
handles.config.C1(3) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function C2_x_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function C2_x_edit_Callback(hObject, eventdata, handles)
handles.config.C2(1) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function C2_y_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function C2_y_edit_Callback(hObject, eventdata, handles)
handles.config.C2(2) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);

% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function C2_z_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function C2_z_edit_Callback(hObject, eventdata, handles)
handles.config.C2(3) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function P1_x_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function P1_x_edit_Callback(hObject, eventdata, handles)
handles.config.P1(1) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function P1_y_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function P1_y_edit_Callback(hObject, eventdata, handles)
handles.config.P1(2) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function P1_z_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function P1_z_edit_Callback(hObject, eventdata, handles)
handles.config.P1(3) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function P2_x_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function P2_x_edit_Callback(hObject, eventdata, handles)
handles.config.P2(1) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function P2_y_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function P2_y_edit_Callback(hObject, eventdata, handles)
handles.config.P2(2) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function P2_z_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function P2_z_edit_Callback(hObject, eventdata, handles)
handles.config.P2(3) = str2num(get(hObject,'String'));
handles = update_topview(handles);
set_axis_limmits(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject, handles);


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function CP_x_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function CP_x_edit_Callback(hObject, eventdata, handles)

x = str2double(get(handles.CP_x_edit,'string'));

if ishandle(handles.topview.current_apex{3})
    handles.config.(handles.topview.current_apex{1})(handles.topview.current_apex{2},1) = x;
    set(handles.topview.current_apex{3},'xdata',x);
    handles = update_topview(handles);
    guidata(hObject,handles);
end


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function CP_y_edit_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function CP_y_edit_Callback(hObject, eventdata, handles)

y = str2double(get(handles.CP_y_edit,'string'));

if ishandle(handles.topview.current_apex{3})
    handles.config.(handles.topview.current_apex{1})(handles.topview.current_apex{2},2) = y;
    set(handles.topview.current_apex{3},'ydata',y);
    handles = update_topview(handles);
    guidata(hObject,handles);
end


% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function Config_popup_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
% --- Executes on selection change in Config_popup.
function Config_popup_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, eventdata, handles)
handles.output = 'Cancel';
uiresume(handles.figure1);


% --------------------------------------------------------------------
% --- Executes on button press in OK_button.
function OK_button_Callback(hObject, eventdata, handles)
handles.output = handles.config;
guidata(hObject,handles);
uiresume(handles.figure1);
%delete(handles.figure1);




% --------------------------------------------------------------------
% --- Executes on button press in Add_button.
function Add_button_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
% --- Executes on button press in Delete_button.
function Delete_button_Callback(hObject, eventdata, handles)

if handles.topview.current_apex{2}~=0
    delete(handles.topview.(handles.topview.current_apex{1})(handles.topview.current_apex{2}));
    handles.topview.(handles.topview.current_apex{1})(handles.topview.current_apex{2}) = [];
    handles.config.(handles.topview.current_apex{1})(handles.topview.current_apex{2},:) = [];
    handles.topview.current_apex{2} = handles.topview.current_apex{2}-1; 
    if handles.topview.current_apex{2} == 0 && isempty(handles.config.(handles.topview.current_apex{1}));
        delete(handles.topview.current_apex{3});
    else
        if handles.topview.current_apex{2} == 0
            handles.topview.current_apex{2} = 1;
        end
        set(handles.topview.current_apex{3},'xdata',...
            handles.config.(handles.topview.current_apex{1})(handles.topview.current_apex{2},1),...
            'ydata',handles.config.(handles.topview.current_apex{1})(handles.topview.current_apex{2},2));
    end
    handles = update_topview(handles);
    guidata(hObject,handles);
end

% --------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)


% --------------------------------------------------------------------
% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)

click_type = get(gcbf,'SelectionType');
switch click_type
    case{'normal'}   % left click
        disp('New Point ... Left click!');
        CurrentPoint = round(mean(get(gca,'CurrentPoint')).*10)./10;
        set(handles.CP_x_edit,'string',num2str(CurrentPoint(1)));
        set(handles.CP_y_edit,'string',num2str(CurrentPoint(2)));
        handles = add_apex(handles,CurrentPoint(1:2),handles.topview.current_apex{1}); 
        set(gcbf,'WindowButtonMotionFcn',{@apex_move_Callback,...
                handles.topview.current_apex{1}});
        set(gcbf,'WindowButtonUpFcn',@ButtonUpFcn);
    case{'extend'}   % Shift - left
    case{'alt'}      % Ctrl - left    
    case{'open'}     % Double click
end

guidata(hObject,handles);


% --------------------------------------------------------------------
% --- Executes on mouse press over apex.
function apex_ButtonDownFcn(hObject, eventdata, wire)%, handles)

handles = guidata(hObject);
click_type = get(gcbf,'SelectionType');
switch click_type
    case{'normal'}   % left click
        disp('Apex ... Left click!');
        
        if ~isfield(handles.config,wire) || isempty(handles.config.(wire))
            if strcmp(wire,'Cwire')
                pos = mean([handles.config.C1(1,:);handles.config.C2(1,:)]);
            else
                pos = mean([handles.config.P1(1,:);handles.config.P2(1,:)]);
            end
            handles.topview.current_apex(1:2) = {wire, 0};
            handles = add_apex(handles,pos(1:2),wire);
            id = 1;
        elseif size(handles.config.(wire),1)==1
            handles.topview.current_apex(1:2) = {wire, 1};
            handles = update_topview(handles);
            id = 1;
        else
            CurrentPoint = mean(get(gca,'CurrentPoint'));
            dist = sqrt((handles.config.(wire)(:,1)-CurrentPoint(1)).^2+...
                (handles.config.(wire)(:,2)-CurrentPoint(2)).^2);
            [temp,id] = min(dist,[],1);
            handles.topview.current_apex(1:2) = {wire, id};
            handles = update_topview(handles);
        end
        if ishandle(handles.topview.current_apex{3})
            set(handles.topview.current_apex{3},'xdata',handles.config.(wire)(id,1),...
                'ydata',handles.config.(wire)(id,2));
        else
            handles.topview.current_apex{3} = line(handles.config.(wire)(id,1),...
                handles.config.(wire)(id,2),'MarkerSize',4, ...
                'markeredgecolor','k','markerfacecolor','k','marker','o',...
                'linewidth',0.5,'ButtonDownFcn',{@apex_ButtonDownFcn, wire});                
        end
            
        %CurrentPoint = round(mean(get(gca,'CurrentPoint')).*10)./10;
        
        %handles = add_apex(handles,CurrentPoint(1:2)); 
        set(gcbf,'WindowButtonMotionFcn',{@apex_move_Callback, wire});
        set(gcbf,'WindowButtonUpFcn',@ButtonUpFcn);
        wire
        set(handles.CP_x_edit,'string',num2str(handles.config.(wire)(id,1)));
        set(handles.CP_y_edit,'string',num2str(handles.config.(wire)(id,2)));
    case{'extend'}   % Shift - left
    case{'alt'}      % Ctrl - left    
    case{'open'}     % Double click
        disp('Apex  Double click!');
end

guidata(hObject,handles);


% --------------------------------------------------------------------
% --- Executes on mouse movement with button pressed on an apex.
function apex_move_Callback(hObject, eventdata, wire)

handles = guidata(hObject);
CurrentPoint = round(mean(get(gca,'CurrentPoint'))*10)/10;
set(handles.CP_x_edit, 'string', num2str(CurrentPoint(1)));
set(handles.CP_y_edit, 'string', num2str(CurrentPoint(2)));

handles.config.(wire)(handles.topview.current_apex{2},1:2) = CurrentPoint(1:2);
if ishandle(handles.topview.current_apex{3})
    set(handles.topview.current_apex{3},'xdata',...
        handles.config.(wire)(handles.topview.current_apex{2},1),...
        'ydata',handles.config.(wire)(handles.topview.current_apex{2},2));
else
    handles.topview.current_apex{3} =...
        line(handles.config.(wire)(handles.topview.current_apex{2},1),...
        handles.config.(wire)(handles.topview.current_apex{2},2),'MarkerSize',4, ...
        'markeredgecolor','k','markerfacecolor','k','marker','o',...
        'linewidth',0.5,'ButtonDownFcn',{@apex_ButtonDownFcn, wire});                
end
set(handles.CP_x_edit,'string',num2str(CurrentPoint(1)));
set(handles.CP_y_edit,'string',num2str(CurrentPoint(2)));


handles = update_topview(handles);

guidata(hObject,handles);


% --------------------------------------------------------------------
% --- Executes on mouse press over electrode.
function electrode_ButtonDownFcn(hObject, eventdata, handles, electrode)

click_type = get(gcbf,'SelectionType');
switch click_type
    case{'normal'}   % left click
        disp(['Electrode ' electrode ' Left click!']);
        set(gcbf,'WindowButtonMotionFcn',{@electrode_move_Callback, electrode});
        set(gcbf,'WindowButtonUpFcn',@ButtonUpFcn);
        set(handles.CP_x_edit, 'string', num2str(handles.config.(electrode)(1)));
        set(handles.CP_y_edit, 'string', num2str(handles.config.(electrode)(2)));
        wire = [electrode(1) 'wire'];
        if ~isempty(handles.config.(wire)) && strcmp(electrode(2),'2')
            apex = size(handles.config.(wire),1);
        else
            apex = 0;
        end
        handles.topview.current_apex(1:2) = {wire, apex};
%        handles = update_topview(handles);
        if ishandle(handles.topview.current_apex{3}) 
            delete(handles.topview.current_apex{3});
            if apex ~= 0
                handles.topview.current_apex{3} = line(handles.config.(wire)(apex,1),...
                   handles.config.(wire)(apex,2),'MarkerSize',4, ...
                   'markeredgecolor','k','markerfacecolor','k','marker','o',...
                   'linewidth',0.5,'ButtonDownFcn',{@apex_ButtonDownFcn, wire});                
           end
        end
    case{'extend'}   % Shift - left
    case{'alt'}      % Ctrl - left    
    case{'open'}     % Double click
        disp(['Electrode ' electrode ' Double click!']);
end
guidata(hObject, handles);

% --------------------------------------------------------------------
% --- Executes on mouse release.
function ButtonUpFcn(hObject, eventdata)

set(gcbf,'WindowButtonMotionFcn','');
set(gcbf,'WindowButtonUpFcn','');

handles = guidata(hObject);
set_axis_limmits(handles);


% --------------------------------------------------------------------
% --- Executes on mouse movement with button pressed on an interface.
function electrode_move_Callback(hObject, eventdata, electrode)

handles = guidata(hObject);
CurrentPoint = mean(get(gca,'CurrentPoint'));
CP = round(CurrentPoint(1:2).*10)./10;
set(handles.CP_x_edit, 'string', num2str(CP(1)));
set(handles.CP_y_edit, 'string', num2str(CP(2)));

handles.config.(electrode)(1:2) = CP;

handles = update_topview(handles);

handles.config.type = 'General Surface Array';
set(handles.Config_popup,'String', handles.config.type);

guidata(hObject,handles);



% --------------------------------------------------------------------
% --- function to update model plot.
function handles = update_topview(handles)

% Update electrode positions
set(handles.C1_x_edit,'string',num2str(handles.config.C1(1)));
set(handles.C1_y_edit,'string',num2str(handles.config.C1(2)));
set(handles.C1_z_edit,'string',num2str(handles.config.C1(3)));
set(handles.topview.C1,'xdata',handles.config.C1(1),...
    'ydata',handles.config.C1(2));

set(handles.C2_x_edit,'string',num2str(handles.config.C2(1)));
set(handles.C2_y_edit,'string',num2str(handles.config.C2(2)));
set(handles.C2_z_edit,'string',num2str(handles.config.C2(3)));
set(handles.topview.C2,'xdata',handles.config.C2(1),...
    'ydata',handles.config.C2(2));

set(handles.P1_x_edit,'string',num2str(handles.config.P1(1)));
set(handles.P1_y_edit,'string',num2str(handles.config.P1(2)));
set(handles.P1_z_edit,'string',num2str(handles.config.P1(3)));
set(handles.topview.P1,'xdata',handles.config.P1(1),...
    'ydata',handles.config.P1(2));

set(handles.P2_x_edit,'string',num2str(handles.config.P2(1)));
set(handles.P2_y_edit,'string',num2str(handles.config.P2(2)));
set(handles.P2_z_edit,'string',num2str(handles.config.P2(3)));
set(handles.topview.P2,'xdata',handles.config.P2(1),...
    'ydata',handles.config.P2(2));

% Update wires
pos1 = handles.config.C1(1:2);
wireunits = length(handles.topview.Cwire);
for k = 1:wireunits
    if isfield(handles.config, 'Cwire') && k<=size(handles.config.Cwire,1)
        pos2 = handles.config.Cwire(k,1:2);
    else
        pos2 = handles.config.C2(1:2);
    end
    set(handles.topview.Cwire(k),'xdata',[pos1(1) pos2(1)],...
        'ydata',[pos1(2) pos2(2)]);
    pos1 = pos2;
end

pos1 = handles.config.P1(1:2);
wireunits = length(handles.topview.Pwire);
for k = 1:wireunits
    if isfield(handles.config, 'Pwire') && k<=size(handles.config.Pwire,1)
        pos2 = handles.config.Pwire(k,1:2);
    else
        pos2 = handles.config.P2(1:2);
    end
    set(handles.topview.Pwire(k),'xdata',[pos1(1) pos2(1)],...
        'ydata',[pos1(2) pos2(2)]);
    pos1 = pos2;
end

h_list = get(handles.axes1,'children');
h_electr = [handles.topview.C1; handles.topview.C2; handles.topview.P1; handles.topview.P2];
h_list = [h_electr; setdiff(h_list, h_electr)];
set(handles.axes1,'children', h_list);

% --------------------------------------------------------------------
% --- function to adjust limmits on axis.

function set_axis_limmits(handles)
minlimmits = min([handles.config.C1;handles.config.C2;handles.config.P1;handles.config.P2],[],1);
maxlimmits = max([handles.config.C1;handles.config.C2;handles.config.P1;handles.config.P2],[],1);

if minlimmits(2)>-20, minlimmits(2) = -20/1.2; end;
if maxlimmits(2)<20, maxlimmits(2) = 20/1.2; end;

set(handles.axes1,'xlim',([minlimmits(1) maxlimmits(1)]+[-0.2 0.2]*diff([minlimmits(1) maxlimmits(1)])));
set(handles.axes1,'ylim',([minlimmits(2) maxlimmits(2)]*1.2));


% --------------------------------------------------------------------
% --- function to add apex to current wire segment.
function handles = add_apex(handles,pos2,wire)

if ~isempty(pos2)
    handles.config.type = 'General Surface Array';
    set(handles.Config_popup,'String', handles.config.type);
        
    %    if strcmp(handles.topview.current_apex{1},'Cwire')
    apex = handles.topview.current_apex{2};
    if isfield(handles.config,wire) && ~isempty(handles.config.(wire))
        if apex == 0
            pos1 = handles.config.([wire(1) '1'])(1,:);
            handles.config.(wire)(2:end+1,:) = ...
               handles.config.(wire)(1:end,:);
        else
            pos1 = handles.config.(wire)(apex,:);
        end
       % handles.config.(wire)(apex+1:end+1,:) = ...
       %     handles.config.(wire)(apex:end,:);
    else
        if strcmp(wire,'Cwire')
            pos1 = handles.config.C1;
        else
            pos1 = handles.config.P1;
        end
    end
    
    handles.topview.(wire)(apex+2:end+1,:) = ...
        handles.topview.(wire)(apex+1:end,:);
    
    handles.config.(wire)(apex+1,1:3) = [pos2 0];
    
    handles.topview.(wire)(apex+1) = line([pos1(1) pos2(1)],...
        [pos1(2) pos2(2)],'MarkerSize',6, ...
        'markeredgecolor','k','markerfacecolor','b','marker','none',...
        'linewidth',2,'color','k', 'ButtonDownFcn',{@apex_ButtonDownFcn, wire});

    
    handles.topview.current_apex(1:2) = {wire, apex+1};
    
    if ishandle(handles.topview.current_apex{3})
        set(handles.topview.current_apex{3},'xdata',handles.config.(wire)(apex+1,1),...
            'ydata',handles.config.(wire)(apex+1,2));
    else
        handles.topview.current_apex{3} = line(handles.config.(wire)(apex+1,1),...
            handles.config.(wire)(apex+1,2),'MarkerSize',4, ...
            'markeredgecolor','k','markerfacecolor','k','marker','o',...
            'linewidth',0.5,'ButtonDownFcn',{@apex_ButtonDownFcn, wire});
    end
        handles = update_topview(handles);
else
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end



