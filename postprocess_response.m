function postprocess_response(d_path)

%This function gets the data of two beads, it then calculates the force
%(the force should be the same on the two beads, and the length change of
%the two beads.
%then we calculate the response by x(w)=res(w)*f(w), so we take teh fourier
%transfor of the force and the length, divide them and select the peak
%value that was already given.
%this gives 2 response function (one for the force on each bead), that
%shoudl be the same.
%This will be done for each data set in the folder, and finally we will get
%the full response function back.

if nargin<1
    d_path=uigetdir('E:\Science\data\elasticity\philippe_no_VASP_with_capping_101029\10');
end

files=dir([d_path,'\*response.mat']);

for j=1:length(files)
    load([d_path,'\',files(j).name]);


    [a1,a2, fr]=get_response_AOD_2beads(data,f,xy_slope,xy_k,cal,rate);
    alpha1(j,:)=a1;
    alpha2(j,:)=a2;
    f_e(j)=f;
end   
    
loglog(f_e,abs(imag(alpha2(:,1))));
hold on
loglog(f_e,abs(imag(alpha1(:,1))),'r');
hold off
