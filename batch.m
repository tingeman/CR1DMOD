function varargout = batch(varargin)

% The batch.m file is related to the batch.fig file, and is used in
% the modeling software CR1Dmod to handle batch calculations of 
% previously defined models.%
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
                   'gui_OpeningFcn', @batch_OpeningFcn, ...
                   'gui_OutputFcn',  @batch_OutputFcn, ...
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


% --- Executes just before batch is made visible.
function batch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to batch (see VARARGIN)

% Choose default command line output for batch
handles.output = hObject;
handles.cancel = 1;

handles.batch = [];
handles.current_model = [];

if nargin > 4
    handles.batch = varargin{2};
    if ~isempty(handles.batch)
        set(handles.Model_listbox, 'string', {handles.batch.name}', 'value', [],...
        'enable', 'on');
    else
        set(handles.Model_listbox, 'string', 'No models added!', 'value', [], ...
            'enable', 'off');
    end
else
    handles.batch = [];
end
if nargin > 3
    handles.current_model = varargin{1};
    if ~isfield(handles.current_model, 'name')
        handles.current_model.name = ['Model ' num2str(length(handles.batch)+1)];
    end
    set(handles.Name_edit, 'String', handles.current_model.name);
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes batch wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = batch_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~handles.cancel
    varargout{1} = handles.current_model;
    varargout{2} = handles.batch;
else
    varargout{1} = [];
    varargout{2} = [];
end

delete(hObject);

% --- Executes during object creation, after setting all properties.
function Name_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Name_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Name_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Name_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Name_edit as text
%        str2double(get(hObject,'String')) returns contents of Name_edit as a double

str = get(hObject,'String');
if ~isempty(str)
    handles.current_model.name = str;
else
    if ~isfield(handles.current_model, 'name') || ...
            isempty(handles.current_model.name)
        handles.current_model.name = ['Model ' num2str(length(handles.batch)+1)];
    end
    set(handles.Name_edit, 'String', handles.current_model.name);
end

guidata(hObject, handles);

% --- Executes on button press in Add_button.
function Add_button_Callback(hObject, eventdata, handles)
% hObject    handle to Add_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.batch)
    handles.batch = handles.current_model;
else
    id = find(ismember(handles.current_model.name, {handles.batch.name}));
    if ~isempty(id)
        button = questdlg('Do you want to replace existing model?',...
            'Model exists!','Yes','No','No');
        if strcmp(button,'Yes')
            handles.batch(id) = handles.current_model;
        elseif strcmp(button,'No')
            return;
        end    
    else
        handles.batch(end+1) = handles.current_model;
    end
end

if ~isempty(handles.batch)
    set(handles.Model_listbox, 'string', {handles.batch.name}',...
        'value', length(handles.batch), 'enable', 'on');
else
    set(handles.Model_listbox, 'string', 'No models added!', 'value', [], ...
        'enable', 'off');
end

guidata(hObject, handles);

% --- Executes on button press in Del_button.
function Del_button_Callback(hObject, eventdata, handles)
% hObject    handle to Del_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

id = get(handles.Model_listbox, 'value');

if ~isempty(id)
    button = questdlg('Delete selected models?',...
        'Delete models?','Yes','No','No');
    if strcmp(button,'Yes')
        handles.batch(id) = [];
    elseif strcmp(button,'No')
        return;
    end    
end
    
if ~isempty(handles.batch)
    set(handles.Model_listbox, 'string', {handles.batch.name}', 'value', [],...
        'enable', 'on');
else
    set(handles.Model_listbox, 'string', 'No models added!', 'value', [], ...
        'enable', 'off');
end

guidata(hObject, handles);


% --- Executes on button press in Up_button.
function Up_button_Callback(hObject, eventdata, handles)
% hObject    handle to Up_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

id = get(handles.Model_listbox, 'value');

if ~isempty(id)
    if id(1) == 1 % disallow selection of first line
        return
        %id = [id(2:end)];
    end

    if ~isempty(id)
        for k = id(:)'
            temp = handles.batch(k-1);
            handles.batch(k-1) = handles.batch(k);
            handles.batch(k) = temp;
        end
        set(handles.Model_listbox, 'string', {handles.batch.name}', 'value', id-1,...
            'enable', 'on');
        
        guidata(hObject, handles);
    end
end


% --- Executes on button press in Down_button.
function Down_button_Callback(hObject, eventdata, handles)
% hObject    handle to Down_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

id = get(handles.Model_listbox, 'value');
if ~isempty(id)
    if id(end) == length(handles.batch) % disallow selection of last line
        return
        %        id = [id(1:end-1)];
    end
    
    if ~isempty(id)
        id = id(end:-1:1);
        for k = id(:)'
            temp = handles.batch(k+1);
            handles.batch(k+1) = handles.batch(k);
            handles.batch(k) = temp;
        end
        id = id(end:-1:1);
        
        set(handles.Model_listbox, 'string', {handles.batch.name}', 'value', id+1,...
            'enable', 'on');
        guidata(hObject, handles);        
    end
end

% --- Executes on button press in OK_button.
function OK_button_Callback(hObject, eventdata, handles)
% hObject    handle to OK_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cancel = 0;
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function Model_listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Model_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in Model_listbox.
function Model_listbox_Callback(hObject, eventdata, handles)
% hObject    handle to Model_listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Model_listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Model_listbox


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

handles.cancel = 1;
guidata(hObject, handles);
uiresume(handles.figure1);
%delete(hObject);


