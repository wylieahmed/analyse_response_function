function [fo,alphax,alphay,trap_stiff,slopes,act_trap]=process_response_function_folder(path)

%the path is the folder where I can find the response data of a given
%set.
%We will load the data, and get the response function for each of the
%frequencies. 
files=dir([path,'\*deformation_response.mat'])
load(([path,'\',files(1).name]))
load(([path,'\results.mat']))
act_trap=find(Traps(:,3)==255);

for j=1:length(files)
    load(([path,'\',files(j).name]))
    fo(j)=f;
    alphax(j)=a_x(act_trap);
    alphay(j)=a_y(act_trap);
    
end

trap_stiff=xy_k;
slopes=xy_slope;