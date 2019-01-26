function varargout = RMS(varargin)
% RMS M-file for RMS.fig
%      RMS, by itself, creates a new RMS or raises the existing
%      singleton*.
%
%      H = RMS returns the handle to a new RMS or the handle to
%      the existing singleton*.
%
%      RMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RMS.M with the given input arguments.
%
%      RMS('Property','Value',...) creates a new RMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RMS_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RMS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RMS

% Last Modified by GUIDE v2.0 12-Feb-2004 13:29:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RMS_OpeningFcn, ...
                   'gui_OutputFcn',  @RMS_OutputFcn, ...
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


% --- Executes just before RMS is made visible.
function RMS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RMS (see VARARGIN)

% Choose default command line output for RMS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RMS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RMS_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function editNsamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNsamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editNsamples_Callback(hObject, eventdata, handles)
% hObject    handle to editNsamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNsamples as text
%        str2double(get(hObject,'String')) returns contents of editNsamples as a double


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global s_i_g_n_a_l;

N = str2double(get(handles.editNsamples,'String'));

[rms,index] = rmsEst(s_i_g_n_a_l,N);

handle_figure = 1111;
figure(handle_figure);
set(handle_figure,'Name','RMS');
set(handle_figure,'NumberTitle','off');

figure(handle_figure);

plot(index,rms);
hold on; plot(s_i_g_n_a_l,'r');

%exporting filtered signal to workspace
assignin('base','r_m_s',rms);
assignin('base','index',index);



% --------------------------------------------------------------------
function varargout = frame1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.frame1.
disp('frame1 Callback not implemented yet.')