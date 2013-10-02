function [fluct,f,psdxy,psdxy_std]=process_fluctuation_folder_v3(path)

%the path is the folder where I can find the fluctuation data of a given
%set.
%We will load the data, and get the average psd and the full linst of PSDs of the files there. Also
%get the std, the list of the slopes, and the active trap to make sure that
%it was not different from the response function
% in the end we store it all in the fluct structure

%this is for testing:
%path='E:\Science\data\response_function\3.28um_beads_rbc\3.28mu_beads_110725\3\multiple_run_series_1003'

files=dir([path,filesep,'fluctuation_data*.mat'])
scan=dir([path,filesep,'scan*.mat'])
load(([path,filesep,files(1).name]))
try 
    act_trap=find(Traps(:,3)==255);
catch 
    load(([path,filesep,'Parameters.mat']))
    Traps=Parameter.Traps;
    act_trap=find(Traps(:,3)==255);
end
    
for i=1:length(files)
    load(([path,'\',files(i).name]))
   
    fluct(i).path=path;
    fluct(i).file=[files(i).name];
    fluct(i).f=fb(act_trap,2:end);
    fluct(i).px=pxb(act_trap,:);
    fluct(i).py=pyb(act_trap,:);
    f(i,:)=fb(act_trap,2:end);
    px(i,:)=pxb(act_trap,:);
    py(i,:)=pyb(act_trap,:);
%     if length(xy_slopes)==0
%         fluct(i).slopes=[1 1 1 1]
% 
%     else
%         fluct(i).slopes=xy_slopes(act_trap,:);
%     end
%     slopes(i,:)=fluct(i).slopes;   
%     
%     load(([path,'\',scan(i).name]))
%     fluct(i).scan_file=[scan(i).name];
%     x(1,:)=x_scan(1,:);
%     x(2,:)=x_scan(4,:)./x_scan(8,:);
%     y(1,:)=y_scan(2,:);
%     y(2,:)=y_scan(6,:)./y_scan(8,:);
%     fluct(i).x_scan=x;
%     fluct(i).y_scan=y;
    
    fluct(i).Traps=Traps
end
f=mean(f,1);

psdx=mean(px,1);
psdx_std=std(px,0,1);

psdy=mean(py,1);
psdy_std=std(py,0,1);

psdxy(1,:)=psdx; 
psdxy(2,:)=psdy; 

psdxy_std(1,:)=psdx_std; 
psdxy_std(2,:)=psdy_std;



