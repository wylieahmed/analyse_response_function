function [f,psdxy,psdxy_std,slopes,act_trap]=process_fluctuation_folder(path)

%the path is the folder where I can find the fluctuation data of a given
%set.
%We will load the data, and get the average psd of the files there. Also
%get the std, the list of the slopes, and the active trap to make sure that
%it was not different from the response function

files=dir([path,'\fluctuation_data*.mat'])
load(([path,'\',files(1).name]))
act_trap=find(Traps(:,3)==255);

for i=1:length(files)
    load(([path,'\',files(i).name]))

    f(i,:)=fl(act_trap,2:end);
    px(i,:)=pxl(act_trap,:);
    py(i,:)=pyl(act_trap,:);
    if length(xy_slopes)==0
        slopes(i,:)=[1 1 1 1]
    else
        slopes(i,:)=xy_slopes(act_trap,:);
    end
    
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


