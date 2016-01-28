function varargout = Fill_Ramp_Rate_Standard(varargin)
% FILL_RAMP_RATE_NEW M-file for Fill_Ramp_Rate_New.fig
%      FILL_RAMP_RATE_NEW, by itself, creates a new FILL_RAMP_RATE_NEW or raises the existing
%      singleton*. 
%
%      H = FILL_RAMP_RATE_NEW returns the handle to a new FILL_RAMP_RATE_NEW or the handle to
%      the existing singleton*.
%
%      FILL_RAMP_RATE_NEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILL_RAMP_RATE_NEW.M with the given input arguments.
%
%      FILL_RAMP_RATE_NEW('Property','Value',...) creates a new FILL_RAMP_RATE_NEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Fill_Ramp_Rate_New_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Fill_Ramp_Rate_New_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Fill_Ramp_Rate_New

% Last Modified by GUIDE v2.5 23-Aug-2015 13:04:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Fill_Ramp_Rate_Standard_OpeningFcn, ...
                   'gui_OutputFcn',  @Fill_Ramp_Rate_Standard_OutputFcn, ...
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


% --- Executes just before Fill_Ramp_Rate_New is made visible.
function Fill_Ramp_Rate_Standard_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Fill_Ramp_Rate_New (see VARARGIN)

% Choose default command line output for Fill_Ramp_Rate_New
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Fill_Ramp_Rate_New wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Fill_Ramp_Rate_Standard_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkboxCD0C.
function checkboxCD0C_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCD0C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCD0C


% --- Executes on button press in checkboxTopOffFueling.
function checkboxTopOffFueling_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTopOffFueling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTopOffFueling


% --- Executes on button press in radiobuttonComm.
function radiobuttonComm_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonComm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonComm

% --- Executes on button press in radiobuttonNonComm.
function radiobuttonNonComm_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonNonComm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonNonComm


% --- Executes on button press in radiobuttonH70a.
function radiobuttonH70a_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonH70a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonH70a


% --- Executes on button press in radiobuttonH70b.
function radiobuttonH70b_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonH70b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonH70b


% --- Executes on button press in radiobuttonH70c.
function radiobuttonH70c_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonH70c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonH70c


% --- Executes on button press in radiobuttonH35a.
function radiobuttonH35a_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonH35a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonH35a


% --- Executes on button press in radiobuttonH35b.
function radiobuttonH35b_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonH35b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonH35b


% --- Executes on button press in radiobuttonT40.
function radiobuttonT40_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonT40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonT40


% --- Executes on button press in radiobuttonT30.
function radiobuttonT30_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonT30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonT30


% --- Executes on button press in radiobuttonT20.
function radiobuttonT20_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonT20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonT20


% --- Executes on button press in checkboxCD10C.
function checkboxCD10C_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCD10C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCD10C


% --- Executes on button press in pushbuttonStartNew.
function pushbuttonStartNew_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonStartNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles = handles = guihandles;
handles = guihandles;
Plot_Fill_Standard_hystep(handles)
close (Fill_Ramp_Rate_Standard);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

delete(hObject);
