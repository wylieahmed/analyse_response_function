function [response]=process_response_function_folder_v2(path)

%the path is the folder where I can find the response data of a given
%set.
%We will load the data, and get the response function for each of the
%frequencies. 
files=dir([path,'\*deformation_response.mat']);
load(([path,'\',files(1).name]))
load(([path,'\results.mat']))
act_trap=find(Traps(:,3)==255);

for j=1:length(files)
    load(([path,'\',files(j).name]))
    fo(j)=f;
    alphax(j)=a_x(act_trap);
    alphay(j)=a_y(act_trap);
    plot(1/s_eff:1/s_eff:length(data)/s_eff,squeeze(data(act_trap,4,:)))
end
trap_stiff=xy_k;
slopes=xy_slope;

response.path=path;
response.f=fo;
response.alphax=alphax;
response.alphay=alphay;
response.trap_stiff=xy_k;
response.slopes=xy_slope;
response.act_trap=act_trap;
response.Traps=Traps;