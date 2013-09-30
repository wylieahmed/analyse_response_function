function [res,f_out]=process_response_data_from_raw_data(path)

%This function will open the raw data of the response function experiments,
%find the active trap, and process to get the real and the imaginary part. 
%It will take the trap stiffness and the position calibration as found in
%the path. Ourput is the response structure as will be used later in
%analysis.

path='\\tweezer-pc\Data\Denis\2012-10-31\response_function_series_1010'
files=dir([path,'\*response.mat']);
load([path,'\parameters.mat']);
act_trap=find(Parameter.Traps(:,3)==255)
AODmu=0.0265;


for j=1:length(files)
    load([path,'\',files(j).name]);
    x_trap=squeeze(data(act_trap,1,:))/AODmu;
    y_trap=squeeze(data(act_trap,2,:))/AODmu;
    power_trap=squeeze(data(act_trap,3,:));
    bead_rel_x=-1*squeeze(data(act_trap,4,:))./squeeze(data(act_trap,8,:))/xy_slope(act_trap,1);
    bead_rel_y=-1*squeeze(data(act_trap,5,:))./squeeze(data(act_trap,8,:))/xy_slope(act_trap,2);

    abs_bead_x=x_trap+bead_rel_x;
    abs_bead_y=y_trap+bead_rel_y;
    
    Fx=1*xy_k(act_trap,1)*bead_rel_x;
    Fy=1*xy_k(act_trap,2)*bead_rel_y;
    
    p=length(abs_bead_x);
% 
     alpha_x=fft(abs_bead_x)./fft(Fx);
     alpha_y=fft(abs_bead_y)./fft(Fy);
%     alpha_x2=fft(x(2,:))./fft(Fx(1,:));
%     alpha_y2=fft(y(2,:))./fft(Fy(1,:));
% 
     fr=s_eff/p*([0:p/2]);
% 
     [m,p_i]=min(abs(fr-f)) ;
     res(j,:)=[alpha_x(p_i),alpha_y(p_i)];
     f_out(j)=f;

end
loglog(f_out,abs(real(res(:,1))));
hold on
loglog(f_out,imag(res(:,1)),'r');
loglog(f_out,imag(1/res(:,1)))
hold off
