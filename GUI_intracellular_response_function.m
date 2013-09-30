function varargout = GUI_intracellular_response_function(varargin)
% GUI_INTRACELLULAR_RESPONSE_FUNCTION MATLAB code for GUI_intracellular_response_function.fig
%      GUI_INTRACELLULAR_RESPONSE_FUNCTION, by itself, creates a new GUI_INTRACELLULAR_RESPONSE_FUNCTION or raises the existing
%      singleton*.
%
%      H = GUI_INTRACELLULAR_RESPONSE_FUNCTION returns the handle to a new GUI_INTRACELLULAR_RESPONSE_FUNCTION or the handle to
%      the existing singleton*.
%
%      GUI_INTRACELLULAR_RESPONSE_FUNCTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_INTRACELLULAR_RESPONSE_FUNCTION.M with the given input arguments.
%
%      GUI_INTRACELLULAR_RESPONSE_FUNCTION('Property','Value',...) creates a new GUI_INTRACELLULAR_RESPONSE_FUNCTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_intracellular_response_function_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_intracellular_response_function_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the axes_response to help GUI_intracellular_response_function

% Last Modified by GUIDE v2.5 24-Sep-2013 11:14:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_intracellular_response_function_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_intracellular_response_function_OutputFcn, ...
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


% --- Executes just before GUI_intracellular_response_function is made visible.
function GUI_intracellular_response_function_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_intracellular_response_function (see VARARGIN)

% Choose default command line output for GUI_intracellular_response_function
handles.output = hObject;

%set the default path
handles.dir='E:\Science\data\response_function\oocytes\2013-06-10';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_intracellular_response_function wait for user axes_response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_intracellular_response_function_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_particle_diameter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_particle_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_particle_diameter as text
%        str2double(get(hObject,'String')) returns contents of edit_particle_diameter as a double
handles=display_current_data(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_particle_diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_particle_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function load_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_response_prefactor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_response_prefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_response_prefactor as text
%        str2double(get(hObject,'String')) returns contents of edit_response_prefactor as a double
handles=display_current_data(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_response_prefactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_response_prefactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function load_folder_Callback(hObject, eventdata, handles)
% hObject    handle to load_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.dir=uigetdir(handles.dir,'Give me the base directory of the response data')
folders=dir([handles.dir,'\response_function*']);
folders(find([folders.isdir]==0))=[];
handles.folders=folders;
handles.current_folder=1;
set(handles.edit_dataset,'String','1');
handles.current_path=folders(handles.current_folder);
handles=load_current_folder(handles);
handles=display_current_data(handles);
guidata(hObject, handles);

%--------------------
function handles=load_current_folder(handles)
%This function loads the data in the current folder

handles.current_folder=str2num(get(handles.edit_dataset,'String'));
folders=handles.folders;
[f,alphax,alphay,trap_stiff,slopes,act_trap]=process_response_function_folder([handles.dir,'\',folders(handles.current_folder).name])

handles.response_f=f;
handles.response_x=alphax;
handles.response_y=alphay;
handles.trap_stiff=trap_stiff;
handles.slopes=slopes;
handles.act_trap=act_trap;

load([handles.dir,'\',folders(handles.current_folder).name,'\histogram_results.mat']);
handles.psd_f=fl(1,2:end);
handles.psd_x=pxl;
handles.psd_y=pyl;

handles.particle_diameter=str2num(get(handles.edit_particle_diameter,'String'))*1e-6;
handles.response_prefactor=str2num(get(handles.edit_response_prefactor,'String'));

%--------------------
function handles=display_current_data(handles)
%This function displays the data in the current folder

handles.particle_diameter=str2num(get(handles.edit_particle_diameter,'String'))*1e-6;
handles.response_prefactor=str2num(get(handles.edit_response_prefactor,'String'));
response_x=handles.response_x;
psd_x=handles.psd_x(handles.act_trap,:);

axes(handles.axes_response)

loglog(handles.response_f,abs(handles.response_prefactor * imag(response_x)));
hold on
loglog(handles.psd_f,abs(pi.*handles.psd_f./4e-21.*psd_x));
hold off

axes(handles.axes_shear_modulus);
eta=1./(3*pi*handles.particle_diameter*2*pi*handles.response_f.*handles.response_prefactor .* abs(imag(handles.response_x)));
loglog(handles.response_f,eta);
hold on
G=1./(6*pi*handles.particle_diameter*handles.response_prefactor .* ((handles.response_x)));
loglog(handles.response_f,abs(real(G)),'+r');
loglog(handles.response_f,abs(imag(G)),'or');
loglog([100:1000],[100:1000].^.75*.1);
hold off
%create dataset to save
folders=handles.folders;


save_data.particle_diameter=handles.particle_diameter;
save_data.response_prefactor=handles.response_prefactor;
save_data.response_f=handles.response_f;
save_data.corr_active_response=abs(handles.response_prefactor * imag(response_x));
save_data.psd_f=handles.psd_f;
save_data.passive_response=abs(pi.*handles.psd_f./4e-21.*psd_x);
save_data.G=G;
save_data.eta=eta;
save_data.base_path=handles.dir;
save_data.folder_path=folders(handles.current_folder).name;

%and now store this structure in the handles.
handles.save_data=save_data;







function edit_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dataset as text
%        str2double(get(hObject,'String')) returns contents of edit_dataset as a double
handles=load_current_folder(handles);
handles=display_current_data(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_dataset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_current.
function save_current_Callback(hObject, eventdata, handles)
% hObject    handle to save_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Here I will save the currently discplayed data in a file that is in the
%current path.
save_data=handles.save_data;
folders=handles.folders;
save([handles.dir,'\',folders(handles.current_folder).name,'\corrected_response.mat'],'save_data');

