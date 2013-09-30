function varargout = inspect_response_function_and_fluctuation_v2(varargin)
% INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION_V2 MATLAB code for inspect_response_function_and_fluctuation_v2.fig
%      INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION_V2, by itself, creates a new INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION_V2 or raises the existing
%      singleton*.
%
%      H = INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION_V2 returns the handle to a new INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION_V2 or the handle to
%      the existing singleton*.
%
%      INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION_V2.M with the given input arguments.
%
%      INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION_V2('Property','Value',...) creates a new INSPECT_RESPONSE_FUNCTION_AND_FLUCTUATION_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inspect_response_function_and_fluctuation_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inspect_response_function_and_fluctuation_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inspect_response_function_and_fluctuation_v2

% Last Modified by GUIDE v2.5 07-Dec-2011 17:50:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inspect_response_function_and_fluctuation_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @inspect_response_function_and_fluctuation_v2_OutputFcn, ...
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


% --- Executes just before inspect_response_function_and_fluctuation_v2 is made visible.
function inspect_response_function_and_fluctuation_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inspect_response_function_and_fluctuation_v2 (see VARARGIN)

% Choose default command line output for inspect_response_function_and_fluctuation_v2
handles.output = hObject;

handles.path='E:\Science\data\response_function\3.28um_beads_rbc\3.28mu_beads_110726\3'
handles.save_path='E:\Science\data\response_function\';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inspect_response_function_and_fluctuation_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inspect_response_function_and_fluctuation_v2_OutputFcn(hObject, eventdata, handles) 
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
xy=get(handles.listbox_xy,'Value');
fluctuation=handles.fluctuation;
slopes=fluctuation(get(handles.listbox_fluct,'Value')).slopes;
set(handles.text_slope,'String',['Actual slope: ',num2str(mean(abs(slopes(:,xy)))),'+-',num2str(std(slopes(:,xy)))]);
set(handles.edit_slope,'String',num2str(mean(abs(slopes(:,xy)))));
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
%handles=update_res(hObject, eventdata, handles)
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
handles=update_fluct(hObject, eventdata, handles)
%handles=update_res(hObject, eventdata, handles)
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

%here I will upen a gui to select the folder with the base data. 
% I will open each of the datasets and preprocess it, and store all in a
% structure which can be later stored.
path=handles.path

path=uigetdir(path);

fluct_files=dir([path,'\multiple_run_series*'])
%make sure I only have directories
fluct_files(find([fluct_files.isdir]~=1))=[];

res_files=dir([path,'\response_function_series*'])
%make sure I only have directories
res_files(find([res_files.isdir]~=1))=[];

%now I write it in the lists. First I build the cell to be shown
%then I will load the full data of interest for each dataset, and store all
%in the fluctuation structure
list_string=[];
clear fluctuation;


for i=1:length(fluct_files)
    fluct_files(i).name
    try
        [fluctuation(i).data,fm,psdxy,psdxy_std,slopes]=process_fluctuation_folder_old_fluct_recording([path,'\',fluct_files(i).name]);
    catch
        [fluctuation(i).data,fm,psdxy,psdxy_std,slopes]=process_fluctuation_folder_new_fluct_recording([path,'\',fluct_files(i).name]);
    end
    fluctuation(i).f=fm;
    fluctuation(i).psdxy=psdxy;
    fluctuation(i).psdxy_std=psdxy_std;
    fluctuation(i).slopes=slopes;  
    fluctuation(i).act_psd=1; 
    fluctuation(i).folder=fluct_files(i).name
    fluctuation(i).path=path
    
end
list_string={fluctuation.folder};
set(handles.listbox_fluct,'String',list_string,'Value',1)
handles.fluctuation=fluctuation;
handles.fluct_files=fluct_files

%Then the response
list_string=[];
clear response;


t=0
for i=1:length(res_files)
    try
        res_int=process_response_function_folder_v2([path,'\',res_files(i).name]);
        res_int.folder=res_files(i).name;
        response(i-t)=res_int;
       
    catch
        t=t+1
    end
end
list_string={response.folder}
set(handles.listbox_res,'String',list_string,'Value',1)
handles.response=response;
handles.res_files=res_files

handles.path=path;
set(handles.text_path,'String',['PATH: ',path])

handles=update_fluct(hObject, eventdata, handles)
%handles=update_res(hObject, eventdata, handles)
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

%[f,psdxy,psdxy_std,slopes,act_trap]=process_fluctuation_folder([handles.path,'\',file]);
fluctuation=handles.fluctuation;
fluct_data=fluctuation(get(handles.listbox_fluct,'Value')).data;

for i=1:length(fluct_data)
    f(i,:)=fluct_data(i).f;
    px(i,:)=fluct_data(i).px;
    py(i,:)=fluct_data(i).py;
    slopes(i,:)=fluct_data(i).slopes;
end

if get(handles.listbox_slope,'Value')==2
    sl=str2num(get(handles.edit_slope,'String'))
    for i=1:size(px,1)        
        px_(i,:)=px(i,:)*(slopes(i,1)/sl)^2;
        py_(i,:)=py(i,:)*(slopes(i,2)/sl)^2;
    end
else
    px_=px;
    py_=py;
end
%(data(i).slopes(1)/sl(i,1))^2


fm=mean(f,1);
psdx=mean(px_,1);
psdx_std=std(px_,0,1); 
psdy=mean(py_,1);
psdy_std=std(py_,0,1);
psdxy(1,:)=psdx; 
psdxy(2,:)=psdy; 
psdxy_std(1,:)=psdx_std; 
psdxy_std(2,:)=psdy_std; 

fluctuation(get(handles.listbox_fluct,'Value')).slopes=mean(abs(slopes))
fluctuation(get(handles.listbox_fluct,'Value')).f=fm;
fluctuation(get(handles.listbox_fluct,'Value')).psdxy=psdxy;
fluctuation(get(handles.listbox_fluct,'Value')).psdxy_std=psdxy_std;
handles.fluctuation=fluctuation;

%now I will read a image to be displayed
handles.image=imread([handles.path,'\',file,'\calibration_end_image.png'])

%now I will display the actual slope
xy=get(handles.listbox_xy,'Value');
set(handles.text_slope,'String',['Actual slope: ',num2str(mean(abs(slopes(:,xy)))),'+-',num2str(std(slopes(:,xy)))]);
%set(handles.edit_slope,'String',num2str(mean(slopes(:,xy))));


% Update handles structure
guidata(hObject, handles);

% % --------------------------------------------------------------------
% function handles=update_res(hObject, eventdata, handles)
% % hObject    handle to path_menu (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% %here I will upen a gui to select the folder with the base data. The dirs I
% %find will be stored in the lists, and can be later selected.
% list=get(handles.listbox_res,'String')
% file=list{get(handles.listbox_res,'Value')};
% 
% [fr,alphax,alphay,trap_stiff,slopes,act_trap]=process_response_function_folder([handles.path,'\',file])
% 
% 
% res_data.fr=fr;
% res_data.alphax=alphax;
% res_data.alphay=alphay;
% res_data.trap_stiff=trap_stiff;
% res_data.slopes=slopes;
% res_data.act_trap=act_trap;
% 
% 
% 
% handles.res_data=res_data;
% 
% %now I will read a image to be displayed
% handles.image=imread([handles.path,'\',file,'\calibration_end_image.png'])
% 
% % Update handles structure
% guidata(hObject, handles);

% --------------------------------------------------------------------
function update_plot(hObject, eventdata, handles)
% hObject    handle to path_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%here I will upen a gui to select the folder with the base data. The dirs I
%find will be stored in the lists, and can be later selected.

%Here I will just read the data from the handles, check if I want to read x
%or y and plot it
r=handles.response;
f=handles.fluctuation;
res_data=r(get(handles.listbox_res,'Value'))
fluct_data=f(get(handles.listbox_fluct,'Value'))
fluctd=fluct_data.data;

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
loglog(res_data.f,abs(imag(res)),'r+');
xlim([.1 1000])
xlabel('Frequency [Hz]')
ylabel('Response [m/N]')

axes(handles.axes_image)
imshow(handles.image)

axes(handles.axes_PSD)
hold off
data=fluct_data.data;
if get(handles.listbox_slope,'Value')==1
    for i=1:length(data)
        sl(i,:)=data(i).slopes;
    end
else
    sl=ones(length(data),2)*str2num(get(handles.edit_slope,'String'))
end

for i=1:length(data)
    if xy==1
        loglog(data(i).f,data(i).px*(data(i).slopes(1)/sl(i,1))^2,'ButtonDownFcn', {@update_active_psd,handles},'Tag', num2str(i))
    else
        loglog(data(i).f,data(i).py*(data(i).slopes(2)/sl(i,2))^2,'ButtonDownFcn', {@update_active_psd,handles},'Tag', num2str(i))
    end
    hold on
end
%now I plot the actual active psd red
    if xy==1
        loglog(data(fluct_data.act_psd).f,data(fluct_data.act_psd).px*(data(fluct_data.act_psd).slopes(1)/sl(fluct_data.act_psd,1))^2,'r')
    else
        loglog(data(fluct_data.act_psd).f,data(fluct_data.act_psd).py*(data(fluct_data.act_psd).slopes(2)/sl(fluct_data.act_psd,2))^2,'r')
    end

axes(handles.axes_slope)
slopes=data(fluct_data.act_psd).slopes
    if xy==1
        scan=data(fluct_data.act_psd).x_scan;
        sl=slopes(1);
    else
        scan=data(fluct_data.act_psd).y_scan;
        sl=slopes(2);
    end
hold off
plot(scan(1,:)-mean(scan(1,:)),scan(2,:)-mean(scan(2,:)));
hold on
plot(scan(1,:)-mean(scan(1,:)),(scan(1,:)-mean(scan(1,:)))*sl,'r');

function edit_slope_Callback(hObject, eventdata, handles)
% hObject    handle to edit_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_slope as text
%        str2double(get(hObject,'String')) returns contents of edit_slope as a double

%If the user change the slope, then we will correct the psd and display it
%in black
% s_new=str2num(get(handles.edit_slope,'String'));
% %the actual slope is:
% fluct_data=handles.fluctuation;
% xy=get(handles.listbox_xy,'Value');
% s_old=(mean(fluct_data.slopes(:,xy)));
handles=update_fluct(hObject, eventdata, handles)
update_plot(hObject, eventdata, handles)
% axes(handles.axes_res)
% hold on
% loglog(fluct_data.f,fluct_data.psdxy(xy,:)*pi.*fluct_data.f./4.1e-21.*(s_old/s_new).^2,'k');


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
%This will remove the actual dataset from the fluctuation
f=handles.fluctuation;
if length(f)>1
    f(get(handles.listbox_fluct,'Value'))=[];
end
handles.fluctuation=f;

list_string={f.folder};
set(handles.listbox_fluct,'String',list_string,'Value',1)

guidata(hObject, handles);
update_plot(hObject, eventdata, handles)

% --- Executes on button press in del_res.
function del_res_Callback(hObject, eventdata, handles)
% hObject    handle to del_res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%This will remove the actual dataset from the fluctuation
r=handles.response;
if length(r)>1
    r(get(handles.listbox_res,'Value'))=[];
end
handles.response=r;

list_string={r.folder};
set(handles.listbox_res,'String',list_string,'Value',1)

guidata(hObject, handles);
update_plot(hObject, eventdata, handles)

% --- Executes on button press in del_psd.
function del_psd_Callback(hObject, eventdata, handles)
% hObject    handle to del_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


fluctuation=handles.fluctuation;
fluct=fluctuation(get(handles.listbox_fluct,'Value'))
fluct_data=fluct.data;
fluct_data(fluct.act_psd)=[]
fluct.act_psd=1
fluct.data=fluct_data;
fluctuation(get(handles.listbox_fluct,'Value'))=fluct
handles.fluctuation=fluctuation;

handles=update_fluct(hObject, eventdata, handles)

guidata(hObject, handles);
update_plot(hObject, eventdata, handles)

% --------------------------------------------------------------------
function save_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to save_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%first we will ask the user for the file location and name
a=strfind(handles.path,'\');

filename=handles.path(a(end-1)+1:end);
filename(strfind(filename,'\'))='-';
[file,path] = uiputfile([handles.save_path,'\',filename,'.mat'],'Save file name');
%then we will save the fluctuation and the response handeles. That should
%do the job.
handles.save_path=path
fluctuation=handles.fluctuation;
response=handles.response;
xy=get(handles.listbox_xy,'Value');
override=get(handles.listbox_slope,'Value');
o_slope=str2num(get(handles.edit_slope,'String'));
save([path,'\',file],'fluctuation','response','xy','override','o_slope');

guidata(hObject, handles);
update_plot(hObject, eventdata, handles)



function update_active_psd(hObject,eventdata,handles)

f=handles.fluctuation;
fluct_data=f(get(handles.listbox_fluct,'Value'))

% get the tag of the selected line
line_hdl = gcbo;
line_tag = get(line_hdl, 'Tag');
f(get(handles.listbox_fluct,'Value')).act_psd=str2num(line_tag);
handles.fluctuation=f;

guidata(hObject, handles);
update_plot(hObject, eventdata, handles)


% --- Executes on selection change in listbox_slope.
function listbox_slope_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_slope contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_slope
handles=update_fluct(hObject, eventdata, handles)
%handles=update_res(hObject, eventdata, handles)
update_plot(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listbox_slope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%first we will ask the user for the file location and name
[file,path] = uigetfile([handles.save_path,'\*.mat'],'Dataset to load');
%then we will save the fluctuation and the response handeles. That should
%do the job.
load([path,'\',file]);

list_string={fluctuation.folder};
set(handles.listbox_fluct,'String',list_string,'Value',1)
handles.fluctuation=fluctuation;

list_string={response.folder}
set(handles.listbox_res,'String',list_string,'Value',1)
handles.response=response;

set(handles.listbox_xy,'Value',xy);
set(handles.listbox_slope,'Value',override);
set(handles.edit_slope,'String',num2str(o_slope));

guidata(hObject, handles);
update_plot(hObject, eventdata, handles)