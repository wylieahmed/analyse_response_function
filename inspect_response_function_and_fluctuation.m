function varargout = inspect_response_function_and_fluctuation(varargin)
% INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION MATLAB code for inspect_response_function_and_fluctuation.fig
%      INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION, by itself, creates a new INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION or raises the existing
%      singleton*.
%
%      H = INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION returns the handle to a new INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION or the handle to
%      the existing singleton*.
%
%      INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION.M with the given input arguments.
%
%      INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION('Property','Value',...) creates a new INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inspect_response_function_and_fluctuation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inspect_response_function_and_fluctuation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inspect_response_function_and_fluctuation

% Last Modified by GUIDE v2.5 29-Nov-2011 14:17:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inspect_response_function_and_fluctuation_OpeningFcn, ...
                   'gui_OutputFcn',  @inspect_response_function_and_fluctuation_OutputFcn, ...
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


% --- Executes just before inspect_response_function_and_fluctuation is made visible.
function inspect_response_function_and_fluctuation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inspect_response_function_and_fluctuation (see VARARGIN)

% Choose default command line output for inspect_response_function_and_fluctuation
handles.output = hObject;

handles.path='E:\Science\data\response_function\3.28um_beads_rbc\3.28mu_beads_110726\3'

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inspect_response_function_and_fluctuation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inspect_response_function_and_fluctuation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_fluct.
function listbox_fluct_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_fluct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_fluct contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_fluct

handles=update_fluct(hObject, eventdata, handles)
update_plot(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listbox_fluct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_fluct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_res.
function listbox_res_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_res contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_res
handles=update_res(hObject, eventdata, handles)
update_plot(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listbox_res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_xy.
function listbox_xy_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_xy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_xy contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_xy
update_plot(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listbox_xy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_xy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function files_menu_Callback(hObject, eventdata, handles)
% hObject    handle to files_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function path_menu_Callback(hObject, eventdata, handles)
% hObject    handle to path_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%here I will upen a gui to select the folder with the base data. The dirs I
%find will be stored in the lists, and can be later selected.
path=handles.path

path=uigetdir(path);

fluct_files=dir([path,'\multiple_run_series*'])
%make sure I only have directories
fluct_files(find([fluct_files.isdir]~=1))=[];

res_files=dir([path,'\response_function_series*'])
%make sure I only have directories
res_files(find([res_files.isdir]~=1))=[];

%now I write it in the lists. First I build the cell to be shown
%First the fluctuation
list_string=[];
for i=1:length(fluct_files)
    curr=fluct_files(i).name;
    list_string{i}=curr;
end
set(handles.listbox_fluct,'String',list_string,'Value',1)
handles.fluct_files=fluct_files

%Then the response
list_string=[];
for i=1:length(res_files)
    curr=res_files(i).name;
    list_string{i}=curr;
end
set(handles.listbox_res,'String',list_string,'Value',1)
handles.res_files=res_files

handles.path=path;
set(handles.text_path,'String',['PATH: ',path])

handles=update_fluct(hObject, eventdata, handles)
handles=update_res(hObject, eventdata, handles)
update_plot(hObject, eventdata, handles)




% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function handles=update_fluct(hObject, eventdata, handles)
% hObject    handle to path_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%here I will upen a gui to select the folder with the base data. The dirs I
%find will be stored in the lists, and can be later selected.

%first I will read the listbox, then I will load the data, and it will be
%stored in a structure that will be added to the handle

list=get(handles.listbox_fluct,'String')
file=list{get(handles.listbox_fluct,'Value')};

[f,psdxy,psdxy_std,slopes,act_trap]=process_fluctuation_folder([handles.path,'\',file]);

fluct_data.f=f;
fluct_data.psdxy=psdxy;
fluct_data.psdxy_std=psdxy_std;
fluct_data.slopes=slopes;
fluct_data.act_trap=act_trap;
handles.fluct_data=fluct_data;

%now I will read a image to be displayed
handles.image=imread([handles.path,'\',file,'\calibration_end_image.png'])

%now I will display the actual slope
xy=get(handles.listbox_xy,'Value');
set(handles.text_slope,'String',['Actual slope: ',num2str(mean(slopes(:,xy))),'+-',num2str(std(slopes(:,xy)))]);
set(handles.edit_slope,'String',num2str(mean(slopes(:,xy))));


% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function handles=update_res(hObject, eventdata, handles)
% hObject    handle to path_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%here I will upen a gui to select the folder with the base data. The dirs I
%find will be stored in the lists, and can be later selected.
list=get(handles.listbox_res,'String')
file=list{get(handles.listbox_res,'Value')};

[fr,alphax,alphay,trap_stiff,slopes,act_trap]=process_response_function_folder([handles.path,'\',file])


res_data.fr=fr;
res_data.alphax=alphax;
res_data.alphay=alphay;
res_data.trap_stiff=trap_stiff;
res_data.slopes=slopes;
res_data.act_trap=act_trap;



handles.res_data=res_data;

%now I will read a image to be displayed
handles.image=imread([handles.path,'\',file,'\calibration_end_image.png'])

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function update_plot(hObject, eventdata, handles)
% hObject    handle to path_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%here I will upen a gui to select the folder with the base data. The dirs I
%find will be stored in the lists, and can be later selected.

%Here I will just read the data from the handles, check if I want to read x
%or y and plot it
res_data=handles.res_data;
fluct_data=handles.fluct_data;
xy=get(handles.listbox_xy,'Value')
if xy==1
    res=res_data.alphax;
else
    res=res_data.alphay;
end

axes(handles.axes_res)
hold off
loglog(fluct_data.f,fluct_data.psdxy(xy,:)*pi.*fluct_data.f/4.1e-21);
hold on
loglog(res_data.fr,abs(imag(res)),'r+');
xlim([.1 1000])
xlabel('Frequency [Hz]')
ylabel('Response [m/N]')

axes(handles.axes_image)
imshow(handles.image)



function edit_slope_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slope as text
%        str2double(get(hObject,'String')) returns contents of edit_slope as a double

%If the user change the slope, then we will correct the psd and display it
%in black
s_new=str2num(get(handles.edit_slope,'String'));
%the actual slope is:
fluct_data=handles.fluct_data;
xy=get(handles.listbox_xy,'Value');
s_old=(mean(fluct_data.slopes(:,xy)));
update_plot(hObject, eventdata, handles)
axes(handles.axes_res)
hold on
loglog(fluct_data.f,fluct_data.psdxy(xy,:)*pi.*fluct_data.f./4.1e-21.*(s_old/s_new).^2,'k');


% --- Executes during object creation, after setting all properties.
function edit_slope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in del_fluct.
function del_fluct_Callback(hObject, eventdata, handles)
% hObject    handle to del_fluct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in del_res.
function del_res_Callback(hObject, eventdata, handles)
% hObject    handle to del_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in del_psd.
function del_psd_Callback(hObject, eventdata, handles)
% hObject    handle to del_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to save_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
