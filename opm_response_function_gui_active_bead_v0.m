function varargout = opm_response_function_gui_active_bead_v0(varargin)
% OPM_RESPONSE_FUNCTION_GUI_ACTIVE_BEAD_V0 M-file for opm_response_function_gui_active_bead_v0.fig
%      OPM_RESPONSE_FUNCTION_GUI_ACTIVE_BEAD_V0, by itself, creates a new OPM_RESPONSE_FUNCTION_GUI_ACTIVE_BEAD_V0 or raises the existing
%      singleton*.
%
%      H = OPM_RESPONSE_FUNCTION_GUI_ACTIVE_BEAD_V0 returns the handle to a new OPM_RESPONSE_FUNCTION_GUI_ACTIVE_BEAD_V0 or the handle to
%      the existing singleton*.
%
%      OPM_RESPONSE_FUNCTION_GUI_ACTIVE_BEAD_V0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPM_RESPONSE_FUNCTION_GUI_ACTIVE_BEAD_V0.M with the given input arguments.
%
%      OPM_RESPONSE_FUNCTION_GUI_ACTIVE_BEAD_V0('Property','Value',...) creates a new OPM_RESPONSE_FUNCTION_GUI_ACTIVE_BEAD_V0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before opm_response_function_gui_active_bead_v0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to opm_response_function_gui_active_bead_v0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help opm_response_function_gui_active_bead_v0

% Last Modified by GUIDE v2.5 20-Jun-2011 12:38:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @opm_response_function_gui_active_bead_v0_OpeningFcn, ...
                   'gui_OutputFcn',  @opm_response_function_gui_active_bead_v0_OutputFcn, ...
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


% --- Executes just before opm_response_function_gui_active_bead_v0 is made visible.
function opm_response_function_gui_active_bead_v0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to opm_response_function_gui_active_bead_v0 (see VARARGIN)

% Choose default command line output for opm_response_function_gui_active_bead_v0
handles.output = hObject;

handles.ddir='E:\Science\data\response_function\3.28mu_beads_110617\3\response_function_series_1003';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes opm_response_function_gui_active_bead_v0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = opm_response_function_gui_active_bead_v0_OutputFcn(hObject, eventdata, handles) 
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
handles=passive_microrheology_one_particle(handles)

   
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
ax=get(handles.listbox1,'value')

%here I will read the passive fluctuation raw data
[anum,param,data_r] = lbv_bin_read_with_parameters([handles.ddir,filesep,'histogram_raw.lvb']);
handles.data_r=data_r;

%now I will update the fields and the plots with the data I have so far...
handles.xy_k=xy_k;
handles.xy_slopes=xy_slope;
set(handles.tab_slopes,'Data',handles.xy_slopes);
set(handles.tab_k,'Data',handles.xy_k);
pos=find(in(:,3)==255)
l=size(in,1);

for i=1:length(res_data)
    f_(i)=res_data(i).f;
    a_x(i)=res_data(i).a_x(pos);
    a_y(i)=res_data(i).a_y(pos);
    j=i;
         %first I want to calculate the displacement. For this I need to calculate
     %the position of each bead in µm anf the substract them
     %first bead
     clear x y Fx Fy dx dy alpha_x1 alpha_y1 alpha_x2 alpha_y2
     data_in=res_data(j).data;
     s_eff=res_data(j).s_eff;
     handles.s_eff=s_eff;
     f=handles.res_data(j).f;
     f_(j)=f;
     for k=1:l
       data=squeeze(data_in(k,:,:));
       x(k,:)=data(1,:)/cal(5)*1e-6-1./(xy_slope(k,1).*1e6).*data(4,:)./data(8,:);
       y(k,:)=data(2,:)/cal(6)*1e-6-1./(xy_slope(k,2).*1e6).*data(6,:)./data(8,:);
       Fx(k,:)=xy_k(k,1)./(xy_slope(1)*1e6).*data(4,:)./data(8,:);
       Fy(k,:)=xy_k(k,2)./(xy_slope(2)*1e6).*data(6,:)./data(8,:);
       
     end
 
     p=length(x);
% 
     alpha_x=fft((x(pos,:)-mean(x,1)))./fft(Fx(pos,:));
     alpha_y=fft((y(pos,:)-mean(y,1)))./fft(Fy(pos,:));
%     alpha_x2=fft(x(2,:))./fft(Fx(1,:));
%     alpha_y2=fft(y(2,:))./fft(Fy(1,:));
% 
     fr=s_eff/p*([0:p/2]);
% 
     [m,p_i]=min(abs(fr-f)) ;
     a1(j,:)=[alpha_x(p_i),alpha_y(p_i)];
     a2(j,:)=[alpha_x(p_i),alpha_y(p_i)];

end
f=f_
axes(handles.res_real);
if (ax==1)
    loglog(f,abs(real(a_x)));
else
    loglog(f,abs(real(a_y)));
end

xlabel('Frequency [Hz]');
ylabel('Response [m/N]');

axes(handles.res_imag);
if (ax==1)
    loglog(f,abs(imag(a_x)));
else
    loglog(f,abs(imag(a_y)));
end
xlabel('Frequency [Hz]');
ylabel('Response [m/N]');

axes(handles.psd_fdt);
if (ax==1)
    loglog(f,4e-21./(pi*f).*abs(imag(a_x)));
else
    loglog(f,4e-21./(pi*f).*abs(imag(a_y)));
end

xlabel('Frequency [Hz]');
ylabel('PSD [m^2/Hz]');

handles.res_data=res_data;
%handles=recalculate_response(hObject, eventdata, handles)

guidata(hObject, handles);    

% % --------------------------------------------------------------------
% function handles=recalculate_response(hObject, eventdata, handles)
% %Here we will recalculate the response data with the correct theory for the
% %two beads, and using the applied QPD factor and stiffness
% %For this I will need to reanalyze the raw data for each frequency:
% cal=handles.cal;
% xy_slopes=get(handles.tab_slopes,'Data');
% xy_k=get(handles.tab_k,'Data');
% handles.xy_slopes=xy_slopes;
% handles.xy_k=xy_k;
% 
% 
% 
% for j=1:length(handles.res_data)
%         %first I want to calculate the displacement. For this I need to calculate
%     %the position of each bead in µm anf the substract them
%     %first bead
%     clear x y Fx Fy dx dy alpha_x1 alpha_y1 alpha_x2 alpha_y2
%     data_in=handles.res_data(j).data;
%     s_eff=handles.res_data(j).s_eff;
%     handles.s_eff=s_eff;
%     f=handles.res_data(j).f;
%     f_(j)=f;
%     for i=1:2
%       data=squeeze(data_in(i,:,:));
%       x(i,:)=data(1,:)/cal(5)*1e-6-1./(xy_slopes(i,1).*1e6).*data(4,:)./data(8,:);
%       y(i,:)=data(2,:)/cal(6)*1e-6-1./(xy_slopes(i,2).*1e6).*data(6,:)./data(8,:);
%       Fx(i,:)=xy_k(i,1)./(xy_slopes(1)*1e6).*data(4,:)./data(8,:);
%       Fy(i,:)=xy_k(i,2)./(xy_slopes(2)*1e6).*data(6,:)./data(8,:);
%     end
% 
%     %now we get the difference in x and y
%     dx=abs(x(1,:)-x(2,:));
%     dy=abs(y(1,:)-y(2,:));
% 
%     p=length(dx);
% 
%     alpha_x1=fft(x(1,:))./fft(Fx(1,:));
%     alpha_y1=fft(y(1,:))./fft(Fy(1,:));
%     alpha_x2=fft(x(2,:))./fft(Fx(1,:));
%     alpha_y2=fft(y(2,:))./fft(Fy(1,:));
% 
%     fr=s_eff/p*([0:p/2]);
% 
%     [m,p_i]=min(abs(fr-f)) ;
%     a1(j,:)=[alpha_x1(p_i),alpha_y1(p_i)];
%     a2(j,:)=[alpha_x2(p_i),alpha_y2(p_i)];
% end
% 
% handles=passive_microrheology(handles);
% 
% 
% %and now use the passive rheology and the previouslt calculated active to
% %get the real response function via the active rheology, and then we are
% %done ;-)
% %The steps are. 
% %1. Extract the passive data at the frequencies where we
% %have the active data. For this I will log bin the data, and the
% %interpolate
% %2. Simply do the translation from apparent to real response function using
% %the equation from Fred.
% 
% %1. Interpolation
% alph=handles.alph;
% Ax=alph.Ax; Ay=alph.Ay;
% f_p=alph.f;
% for pos=1:2
%     [fl(pos,:),Axl(pos,:)]=log_binning(f_p,Ax(pos,:),200);
%     [fl(pos,:),Ayl(pos,:)]=log_binning(f_p,Ay(pos,:),200);
%     Ax_i(pos,:)=interp1(fl(pos,:),Axl(pos,:),f_);
%     Ay_i(pos,:)=interp1(fl(pos,:),Ayl(pos,:),f_);
% end
% %2 Real active response data
% ax_a=a2(:,1)'./((1-xy_k(1,1)*Ax_i(1,:)).*(1-xy_k(2,1)*Ax_i(2,:))-xy_k(1,1)*xy_k(2,1)*a2(:,1)');
% ay_a=a2(:,2)'./((1-xy_k(1,2)*Ay_i(1,:)).*(1-xy_k(2,2)*Ay_i(2,:))-xy_k(1,2)*xy_k(2,2)*a2(:,2)');
% 
% alph.ax_a=ax_a;
% alph.ay_a=ay_a;
% 
% %Now we put the active one and two particel response in the alph structure
% alph.amr_tp.ax_a=ax_a;
% alph.amr_tp.ay_a=ay_a;
% alph.amr_tp.f=f_;
% 
% alph.amr_op.a1x_a=a1(:,1)./(1-xy_k(1,1)*a1(:,1));
% alph.amr_op.a1y_a=a1(:,2)./(1-xy_k(1,2)*a1(:,2));
% alph.amr_op.f=f_;
% 
% 
% handles.alph=alph;
% 
% 
% 
% axes(handles.res_real);
% %loglog(f_,abs(real(alph.amr)),f_,abs(real(a2)));
% loglog(alph.amr_tp.f,abs(real(-alph.amr_tp.ax_a)),alph.pmr_tp.f,abs(real(-alph.pmr_tp.ax_p)));
% xlabel('Frequency [Hz]');
% ylabel('Response [m/N]');
% 
% axes(handles.res_imag);
% %loglog(f_,abs(imag(a1)),f_,abs(imag(a2)));
% loglog(alph.amr_tp.f,abs(imag(-alph.amr_tp.ax_a)),alph.pmr_tp.f,abs(imag(-alph.pmr_tp.ax_p)));
% xlabel('Frequency [Hz]');
% ylabel('Response [m/N]');
% 
% axes(handles.psd_fdt);
% %loglog(f_,4e-20./(pi*[f_;f_]').*abs(imag(a1)));%,f_,4e-20./(pi*f_).*abs(imag(a1(2,:)')));
% loglog(alph.amr_op.f,abs(imag(-alph.amr_op.a1x_a)),alph.pmr_op.f,abs(imag(-alph.pmr_op.a1x_p)));
% 
% xlabel('Frequency [Hz]');
% ylabel('PSD [m^2/Hz]');
    
%---------------------------------------------------------------
function handles=passive_microrheology_one_particle(handles)

% The path is simple:
%1: calculate the PSD of active particle, and of the interparticle distance
%2: Use the FDT to get the imaginary part of the response function
%3: Use the Kramers Kronig theorem to get the real part


%Hense, the input here needs to be the raw fluctuation data, the trap
%stiffness and the calibration factor, the scanrate

corr=str2num(get(handles.corr1,'String'));
data_r=handles.data_r;
rdata=handles.rdata;
o=rdata.Data.Parameters.Root.offtime_in_the_begining.value;
xy_slopes=handles.xy_slopes;
xy_k=handles.xy_k;
ax=get(handles.listbox1,'value')

hsrt=handles.rdata.Data.Parameters.Root.Sampling_rate.value;

srt=hsrt/(handles.rdata.Data.Parameters.Root.Average_per_Trap.value+handles.rdata.Data.Parameters.Root.AOD_settling_time__micros.value*1e-6*hsrt);
srt=srt/size(handles.in,1);

cal=handles.cal;
in=handles.in;
data=data_r(:,size(in,1)*o+1:end);
%1: PSD of active particle 
        pos=find(in(:,3)==255)
        l=size(in,1);
            sh=[pos:l:length(data)];
            xi=data(1,sh)./data(3,sh)/(xy_slopes(pos,1)*1e6);
            x=xi;
            y=data(2,sh)./data(3,sh)/(xy_slopes(pos,2)*1e6);
            [f_p,px]=power_sd(x,srt);
            [f_p,py]=power_sd(y,srt);
            [fl,pxl]=log_binning(f_p,px,200);
            [fl,pyl]=log_binning(f_p,py,200);
       

            
 %2: Use the FDT to get the imaginary part of the response function
 %I will store the apparent response function an A, where the first index
 %is for the first, second an interparticle distance, and the cesond
 %dimension is for the frequencies
 f=f_p;
 Axii=pi*f.*px/4e-21;
 Ayii=pi*f.*py/4e-21;
 
 %3: use the kramers kronig relation to get the real part
 
    Ax=kk_rel(Axii);
    Ay=kk_rel(Ayii);
    
    clear i
    
ax_p=Ax+i*Axii;
ay_p=Ay+i*Ayii;



%finally we will do a log binning for convenience, each dataset will be
%stored as the log binned version. Only Ax and AY will be in the full and
%the log binned version
[fl,ax_pl]=log_binning(f(1,:),ax_p,200);
[fl,ay_pl]=log_binning(f(1,:),ay_p,200);




%And nowI put it all in the handles structure that will be returned
%I will assemble a alpha structure that will late contain all the info
%The final structure will be seperated into pactive and passive, and one or
%2 partle micrhorheology, and apparent response function A
%so. alph.pmr_op (passive_micrhorheology_onepartile), alph.pmr_tp
%(passive_micrhorheology_twoparticle)
alph.ax_p=ax_p;
alph.ay_p=ay_p;
alph.ax_pl=ax_pl;
alph.ay_pl=ay_pl;
alph.f=f(1,:);
alph.fl=fl;

%and now plot the psd
axes(handles.psd_fdt);
hold on
if (ax==1)
    loglog(fl,pxl*corr,'r');
else
    loglog(fl,pyl*corr,'r');
end
hold off
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


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corr1_Callback(hObject, eventdata, handles)
% hObject    handle to corr1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corr1 as text
%        str2double(get(hObject,'String')) returns contents of corr1 as a double


% --- Executes during object creation, after setting all properties.
function corr1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corr1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
