function varargout = response_function_gui_v0(varargin)
% RESPONSE_FUNCTION_GUI_V0 M-file for response_function_gui_v0.fig
%      RESPONSE_FUNCTION_GUI_V0, by itself, creates a new RESPONSE_FUNCTION_GUI_V0 or raises the existing
%      singleton*.
%
%      H = RESPONSE_FUNCTION_GUI_V0 returns the handle to a new RESPONSE_FUNCTION_GUI_V0 or the handle to
%      the existing singleton*.
%
%      RESPONSE_FUNCTION_GUI_V0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESPONSE_FUNCTION_GUI_V0.M with the given input arguments.
%
%      RESPONSE_FUNCTION_GUI_V0('Property','Value',...) creates a new RESPONSE_FUNCTION_GUI_V0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before response_function_gui_v0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to response_function_gui_v0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help response_function_gui_v0

% Last Modified by GUIDE v2.5 14-Jan-2011 13:04:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @response_function_gui_v0_OpeningFcn, ...
                   'gui_OutputFcn',  @response_function_gui_v0_OutputFcn, ...
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


% --- Executes just before response_function_gui_v0 is made visible.
function response_function_gui_v0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to response_function_gui_v0 (see VARARGIN)

% Choose default command line output for response_function_gui_v0
handles.output = hObject;

handles.ddir='E:\Science\data\response_function\1µm_bead_in_water\110113\2\response_function_series_1000';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes response_function_gui_v0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = response_function_gui_v0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function qpd_cal_Callback(hObject, eventdata, handles)
% hObject    handle to qpd_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of qpd_cal as text
%        str2double(get(hObject,'String')) returns contents of qpd_cal as a double
handles=recalculate_response(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function qpd_cal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to qpd_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function trap_stiff_Callback(hObject, eventdata, handles)
% hObject    handle to trap_stiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trap_stiff as text
%        str2double(get(hObject,'String')) returns contents of trap_stiff as a double
handles=recalculate_response(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function trap_stiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trap_stiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ddir=handles.ddir;
handles.ddir=uigetdir(ddir,'Please specify the parent dir');
%First I will read the mat files, which should have all the date required
%to redo the analysis with different calibration values.
handles=load_response_data(hObject, eventdata, handles);


   
handles.ddir=ddir;
guidata(hObject, handles);

% --------------------------------------------------------------------
function handles=load_response_data(hObject, eventdata, handles)
%Here I will oad all the  *deformation_response.mat files
% Which are supposed to contain the info required to reconstruct the
% Response with new calibration values.
files=dir([handles.ddir,filesep,'*deformation_response.mat']);

%now I open each of the files and store the important data in a structure
%callred res_data
for j=1:length(files)
    load([handles.ddir,filesep,files(j).name])
    res_data(j).data=data;
    res_data(j).scan_path=scan_path;
    res_data(j).repead_oss=repead_oss;
    res_data(j).f=f;
    res_data(j).a_x=a_x;
    res_data(j).a_y=a_y;
    res_data(j).s_eff=s_eff;
end


%new I will read the results tdms file
[rdata]=convertTDMS(0,[handles.ddir,filesep,'results.tdms']);

for i=1:length(rdata.Data.MeasuredData)-14
 in(i,:)=rdata.Data.MeasuredData(14+i).Data;
end 

handles.in=in;
handles.cal=cal;
handles.rdata=rdata;

%here I will read the passive fluctuation raw data
[anum,param,data_r] = lbv_bin_read_with_parameters([handles.ddir,filesep,'histogram_raw.lvb']);
handles.data_r=data_r;

%now I will update the fields and the plots with the data I have so far...
handles.xy_k=xy_k;
handles.xy_slopes=xy_slope(1:2,1:2);
set(handles.tab_slopes,'Data',handles.xy_slopes);
set(handles.tab_k,'Data',handles.xy_k);

for i=1:length(res_data)
    f(i)=res_data(i).f;
    a_x(i)=res_data(i).a_x(1);
    a_y(i)=res_data(i).a_y(1);
end
axes(handles.res_real);
loglog(f,abs(real(a_x)),f,abs(real(a_y)));
xlabel('Frequency [Hz]');
ylabel('Response [m/N]');

axes(handles.res_imag);
loglog(f,abs(imag(a_x)),f,abs(imag(a_y)));
xlabel('Frequency [Hz]');
ylabel('Response [m/N]');

axes(handles.psd_fdt);
loglog(f,4e-20./(pi*f).*abs(imag(a_x)),f,4e-20./(pi*f).*abs(imag(a_y)));
xlabel('Frequency [Hz]');
ylabel('PSD [m^2/Hz]');

handles.res_data=res_data;
handles=recalculate_response(hObject, eventdata, handles)

guidata(hObject, handles);    

% --------------------------------------------------------------------
function handles=recalculate_response(hObject, eventdata, handles)
%Here we will recalculate the response data with the correct theory for the
%two beads, and using the applied QPD factor and stiffness
%For this I will need to reanalyze the raw data for each frequency:
cal=handles.cal;
xy_slopes=get(handles.tab_slopes,'Data');
xy_k=get(handles.tab_k,'Data');
handles.xy_slopes=xy_slopes;
handles.xy_k=xy_k;



for j=1:length(handles.res_data)
        %first I want to calculate the displacement. For this I need to calculate
    %the position of each bead in µm anf the substract them
    %first bead
    clear x y Fx Fy dx dy alpha_x1 alpha_y1 alpha_x2 alpha_y2
    data_in=handles.res_data(j).data;
    s_eff=handles.res_data(j).s_eff;
    handles.s_eff=s_eff;
    f=handles.res_data(j).f;
    f_(j)=f;
    for i=1:2
      data=squeeze(data_in(i,:,:));
      x(i,:)=data(1,:)/cal(5)*1e-6-1./(xy_slopes(i,1).*1e6).*data(4,:)./data(8,:);
      y(i,:)=data(2,:)/cal(6)*1e-6-1./(xy_slopes(i,2).*1e6).*data(6,:)./data(8,:);
      Fx(i,:)=xy_k(i,1)./(xy_slopes(1)*1e6).*data(4,:)./data(8,:);
      Fy(i,:)=xy_k(i,2)./(xy_slopes(2)*1e6).*data(6,:)./data(8,:);
    end

    %now we get the difference in x and y
    dx=abs(x(1,:)-x(2,:));
    dy=abs(y(1,:)-y(2,:));

    p=length(dx);

    alpha_x1=fft(x(1,:))./fft(Fx(1,:));
    alpha_y1=fft(y(1,:))./fft(Fy(1,:));
    alpha_x2=fft(x(2,:))./fft(Fx(1,:));
    alpha_y2=fft(y(2,:))./fft(Fy(1,:));

    fr=s_eff/p*([0:p/2]);

    [m,p_i]=min(abs(fr-f)) ;
    a1(j,:)=[alpha_x1(p_i),alpha_y1(p_i)];
    a2(j,:)=[alpha_x2(p_i),alpha_y2(p_i)];
end

handles=passive_microrheology(handles);


%and now use the passive rheology and the previouslt calculated active to
%get the real response function via the active rheology, and then we are
%done ;-)
%The steps are. 
%1. Extract the passive data at the frequencies where we
%have the active data. For this I will log bin the data, and the
%interpolate
%2. Simply do the translation from apparent to real response function using
%the equation from Fred.

%1. Interpolation
alph=handles.alph;
Ax=alph.Ax; Ay=alph.Ay;
f_p=alph.f;
for pos=1:2
    [fl(pos,:),Axl(pos,:)]=log_binning(f_p,Ax(pos,:),200);
    [fl(pos,:),Ayl(pos,:)]=log_binning(f_p,Ay(pos,:),200);
    Ax_i(pos,:)=interp1(fl(pos,:),Axl(pos,:),f_);
    Ay_i(pos,:)=interp1(fl(pos,:),Ayl(pos,:),f_);
end
%2 Real active response data
ax_a=a2(:,1)'./((1-xy_k(1,1)*Ax_i(1,:)).*(1-xy_k(2,1)*Ax_i(2,:))-xy_k(1,1)*xy_k(2,1)*a2(:,1)');
ay_a=a2(:,2)'./((1-xy_k(1,2)*Ay_i(1,:)).*(1-xy_k(2,2)*Ay_i(2,:))-xy_k(1,2)*xy_k(2,2)*a2(:,2)');

alph.ax_a=ax_a;
alph.ay_a=ay_a;

%Now we put the active one and two particel response in the alph structure
alph.amr_tp.ax_a=ax_a;
alph.amr_tp.ay_a=ay_a;
alph.amr_tp.f=f_;

alph.amr_op.a1x_a=a1(:,1)./(1-xy_k(1,1)*a1(:,1));
alph.amr_op.a1y_a=a1(:,2)./(1-xy_k(1,2)*a1(:,2));
alph.amr_op.f=f_;


handles.alph=alph;



axes(handles.res_real);
%loglog(f_,abs(real(alph.amr)),f_,abs(real(a2)));
loglog(alph.amr_tp.f,abs(real(-alph.amr_tp.ax_a)),alph.pmr_tp.f,abs(real(-alph.pmr_tp.ax_p)));
xlabel('Frequency [Hz]');
ylabel('Response [m/N]');

axes(handles.res_imag);
%loglog(f_,abs(imag(a1)),f_,abs(imag(a2)));
loglog(alph.amr_tp.f,abs(imag(-alph.amr_tp.ax_a)),alph.pmr_tp.f,abs(imag(-alph.pmr_tp.ax_p)));
xlabel('Frequency [Hz]');
ylabel('Response [m/N]');

axes(handles.psd_fdt);
%loglog(f_,4e-20./(pi*[f_;f_]').*abs(imag(a1)));%,f_,4e-20./(pi*f_).*abs(imag(a1(2,:)')));
loglog(alph.amr_op.f,abs(imag(-alph.amr_op.a1x_a)),alph.pmr_op.f,abs(imag(-alph.pmr_op.a1x_p)));

xlabel('Frequency [Hz]');
ylabel('PSD [m^2/Hz]');
    
%---------------------------------------------------------------
function handles=passive_microrheology(handles)
%ok, now I need to get from these data the response function. I will follow
%the paper by Mizuno (Macrorheology 2008,41,7194-1202
%This means that I should first get the single particle response function
%from the passive data, This will require each trap stiffness, and the
%fluctuation of each data.
% The path is simple:
%1: calculate the PSD of each particle, and of the interparticle distance
%2: Use the FDT to get the imaginary part of the response function
%3: Use the Kramers Kronig theorem to get the real part
%4: Use the formular from Fred Mackintosh to get the real response function
%5: Additional use the data of each particle to get the one particle
%response function

%Hense, the input here needs to be the raw fluctuation data, the trap
%stiffness and the calibration factor, the scanrate

data_r=handles.data_r;
rdata=handles.rdata;
o=rdata.Data.Parameters.Root.offtime_in_the_begining.value;
xy_slopes=handles.xy_slopes;
xy_k=handles.xy_k;

hsrt=handles.rdata.Data.Parameters.Root.Sampling_rate.value;

srt=hsrt/(handles.rdata.Data.Parameters.Root.Average_per_Trap.value+handles.rdata.Data.Parameters.Root.AOD_settling_time__micros.value*1e-6*hsrt);
srt=srt/size(handles.in,1);

cal=handles.cal;
in=handles.in;
data=data_r(:,size(in,1)*o+1:end);
%1: PSD of each particle and of the interparticle distance

        %this is for the 2 traps
        for pos=1:2
            sh=[pos:2:length(data)];
            xi=data(1,sh)./data(3,sh)/(xy_slopes(pos,1)*1e6);
            x(pos,:)=xi;
            y(pos,:)=data(2,sh)./data(3,sh)/(xy_slopes(pos,2)*1e6);
            [f_p(pos,:),px(pos,:)]=power_sd(x(pos,:),srt);
            [f_p(pos,:),py(pos,:)]=power_sd(y(pos,:),srt);
            [fl(pos,:),pxl(pos,:)]=log_binning(f_p(pos,:),px(pos,:),200);
            [fl(pos,:),pyl(pos,:)]=log_binning(f_p(pos,:),py(pos,:),200);
        end
        %and now I calculate for the interparticle distance
        pos=3;
        x(pos,:)=(in(1,1)-in(2,1))/cal(5)*1e-6-(x(1,:)-x(2,:));
        y(pos,:)=(in(1,2)-in(2,2))/cal(6)*1e-6-(y(1,:)-y(2,:));

        %Then we get the interparticle powerspectrum
            [f_p(pos,:),px(pos,:)]=power_sd(x(pos,:),srt);
            [f_p(pos,:),py(pos,:)]=power_sd(y(pos,:),srt);
            [fl(pos,:),pxl(pos,:)]=log_binning(f_p(pos,:),px(pos,:),200);
            [fl(pos,:),pyl(pos,:)]=log_binning(f_p(pos,:),py(pos,:),200);
            
 %2: Use the FDT to get the imaginary part of the response function
 %I will store the apparent response function an A, where the first index
 %is for the first, second an interparticle distance, and the cesond
 %dimension is for the frequencies
 f=f_p;
 Axii=pi*f.*px/4e-20;
 Ayii=pi*f.*py/4e-20;
 
 %3: use the kramers kronig relation to get the real part
 
for i=1:3
    Ax(i,:)=kk_rel(Axii(i,:));
    Ay(i,:)=kk_rel(Ayii(i,:));
end

%4: Use the formular from Fred Mackintosh to get the real response function
%For this I need the trap stiffness of each trap
ax_p=Ax(3,:)./((1-xy_k(1,1)*Ax(1))*(1-xy_k(2,1)*Ax(2,:))-xy_k(1,1)*xy_k(2,1)*Ax(3,:));
ay_p=Ay(3,:)./((1-xy_k(1,2)*Ay(1))*(1-xy_k(2,2)*Ay(2,:))-xy_k(1,2)*xy_k(2,2)*Ay(3,:));

%5: Additional use the data of each particle to get the one particle
%response function

a1x_p=Ax(1,:)./(1-xy_k(1,1)*Ax(1,:));
a1y_p=Ay(1,:)./(1-xy_k(1,2)*Ay(1,:));

a2x_p=Ax(2,:)./(1-xy_k(2,1)*Ax(2,:));
a2y_p=Ay(2,:)./(1-xy_k(2,2)*Ay(2,:));

%finally we will do a log binning for convenience, each dataset will be
%stored as the log binned version. Only Ax and AY will be in the full and
%the log binned version
[fl,ax_pl]=log_binning(f(1,:),ax_p,200);
[fl,ay_pl]=log_binning(f(1,:),ay_p,200);
[fl,a1x_pl]=log_binning(f(1,:),a1x_p,200);
[fl,a1y_pl]=log_binning(f(1,:),a1y_p,200);
[fl,a2x_pl]=log_binning(f(1,:),a2x_p,200);
[fl,a2y_pl]=log_binning(f(1,:),a2y_p,200);



%And nowI put it all in the handles structure that will be returned
%I will assemble a alpha structure that will late contain all the info
%The final structure will be seperated into pactive and passive, and one or
%2 partle micrhorheology, and apparent response function A
%so. alph.pmr_op (passive_micrhorheology_onepartile), alph.pmr_tp
%(passive_micrhorheology_twoparticle)
alph.Ax=Ax;
alph.Ay=Ay;
alph.pmr_tp.ax_p=ax_pl;
alph.pmr_tp.ay_p=ay_pl;
alph.pmr_tp.f=fl;
alph.pmr_op.a1x_p=a1x_pl;
alph.pmr_op.a1y_p=a1y_pl;
alph.pmr_op.a2x_p=a2x_pl;
alph.pmr_op.a2y_p=a2y_pl;
alph.pmr_op.f=fl;
alph.f=f(1,:);
alph.fl=fl;
handles.alph=alph;

%-------------------------------------------------------------------------
function chi=kk_rel(chi_ii);
%This uses the method of Hilbert transofrmation by Pietro
inter=-1*chi_ii(end:-1:1);
inter=[inter chi_ii];
hil=hilbert(inter);
hil=hil(end/2+1:end);
chi=imag(hil)+i*chi_ii;


% --- Executes when entered data in editable cell(s) in tab_slopes.
function tab_slopes_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tab_slopes (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

handles=recalculate_response(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes when entered data in editable cell(s) in tab_k.
function tab_k_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to tab_k (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

handles=recalculate_response(hObject, eventdata, handles);
guidata(hObject, handles);
