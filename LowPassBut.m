function varargout = LowPassBut(varargin)
% LOWPASSBUT M-file for LowPassBut.fig
%      LOWPASSBUT, by itself, creates a new LOWPASSBUT or raises the existing
%      singleton*.
%
%      H = LOWPASSBUT returns the handle to a new LOWPASSBUT or the handle to
%      the existing singleton*.
%
%      LOWPASSBUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOWPASSBUT.M with the given input arguments.
%
%      LOWPASSBUT('Property','Value',...) creates a new LOWPASSBUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LowPassBut_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LowPassBut_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LowPassBut

% Last Modified by GUIDE v2.5 12-Nov-2003 11:08:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LowPassBut_OpeningFcn, ...
                   'gui_OutputFcn',  @LowPassBut_OutputFcn, ...
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


% --- Executes just before LowPassBut is made visible.
function LowPassBut_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LowPassBut (see VARARGIN)

% Choose default command line output for LowPassBut
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LowPassBut wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LowPassBut_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonProcess.
function pushbuttonProcess_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonProcess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes on button press in pushbuttonApplyFilter.
function pushbuttonApplyFilter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonApplyFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global s_i_g_n_a_l;


fsr = str2double(get(handles.editSampFreq,'String'));
fc = str2double(get(handles.edit_cutoffFreq,'String'));
n = str2double(get(handles.edit_order,'String'));

[filterCoef,ybf,y,outSignal]= Butfilter (s_i_g_n_a_l,fsr,fc,n);

%exporting envelopes to workspace
assignin('base','ButFilteredSig',outSignal);
assignin('base','FilterCoefficients',filterCoef);


% --- Executes during object creation, after setting all properties.
function edit_sampFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sampFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_sampFreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sampFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sampFreq as text
%        str2double(get(hObject,'String')) returns contents of edit_sampFreq as a double


% --- Executes during object creation, after setting all properties.
function editSampFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSampFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editSampFreq_Callback(hObject, eventdata, handles)
% hObject    handle to editSampFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSampFreq as text
%        str2double(get(hObject,'String')) returns contents of editSampFreq as a double


% --- Executes during object creation, after setting all properties.
function edit_cutoffFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cutoffFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_cutoffFreq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cutoffFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cutoffFreq as text
%        str2double(get(hObject,'String')) returns contents of edit_cutoffFreq as a double


% --- Executes during object creation, after setting all properties.
function edit_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit_order_Callback(hObject, eventdata, handles)
% hObject    handle to edit_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_order as text
%        str2double(get(hObject,'String')) returns contents of edit_order as a double


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbuttonApplyFilter.
function pushbuttonApplyFilter_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbuttonApplyFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


