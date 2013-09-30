function [res,f_out,res_all_x,res_all_y]=process_response_data_for_opposing_beads_from_raw_data(path)

%This function will open the raw data of the response function experiments,
%find the active trap, and try to get the opposing bead. Then get the real and the imaginary part of the active part and of the correlation between the other beads. 
%It will take the trap stiffness and the position calibration as found in
%the path. Output is the response structure as will be used later in
%analysis.
%Furthermore it will output a response structure that gives the correlation
%mixes between the different beads with the active bead. The datasets are
%either the absolute value of the active bead, or the distance values with
%respect to that bead

%path='E:\Science\data\response_function\3.28um_beads_rbc\3.28mu_beads_110726\1\response_function_series_1000'
files=dir([path,'\*response.mat']);
load([path,'\parameters.mat']);
act_trap=find(Parameter.Traps(:,3)==255)
n_traps=size(Parameter.Traps,1);
AODmu=0.0265;


for j=1:length(files)
    %load the data
    load([path,'\',files(j).name]);
    %here I extract the trap and the bead positions of each of the traps
    clear x_trap y_trap power_trap bead_rel_x bead_rel_y abs_bead_x abs_bead_y Fx Fy alpha_x alpha_y 
    for i=1:n_traps
        x_trap(i,:)=squeeze(data(i,1,:))/AODmu;
        y_trap(i,:)=squeeze(data(i,2,:))/AODmu;
        power_trap(i,:)=squeeze(data(i,3,:));
        bead_rel_x(i,:)=-1*squeeze(data(i,4,:))./squeeze(data(i,8,:))/xy_slope(i,1);
        bead_rel_y(i,:)=-1*squeeze(data(i,5,:))./squeeze(data(i,8,:))/xy_slope(i,2);
        abs_bead_x(i,:)=x_trap(i,:)+bead_rel_x(i,:);
        abs_bead_y(i,:)=y_trap(i,:)+bead_rel_y(i,:);
        Fx(i,:)=1*xy_k(act_trap,1)*bead_rel_x(i,:);
        Fy(i,:)=1*xy_k(act_trap,2)*bead_rel_y(i,:);    
    end
    p=length(abs_bead_x);
    for i=1:n_traps
        if i==act_trap
            %this is now the case of the active bead       
            alpha_x(i,:)=fft(abs_bead_x(act_trap,:))./fft(Fx(act_trap,:));
            alpha_y(i,:)=fft(abs_bead_y(act_trap,:))./fft(Fy(act_trap,:));
        else
            %here I need to get the correlation between the current and the
            %active bead. 
            %in fact I will calculate the distance (x,y) between the two
            %beads. Then I will take the fft of the distance, and the fftof
            %the force applied on the active bead. This way I might
            %introduce an artefact as there are two other beads that should
            %be taken into account. 
            alpha_x(i,:)=fft(abs_bead_x(act_trap,:)-abs_bead_x(i,:))./fft(Fx(act_trap,:));
            alpha_y(i,:)=fft(abs_bead_y(act_trap,:)-abs_bead_y(i,:))./fft(Fy(act_trap,:));
        end
    end
        

%     alpha_x2=fft(x(2,:))./fft(Fx(1,:));
%     alpha_y2=fft(y(2,:))./fft(Fy(1,:));
% 
     fr=s_eff/p*([0:p/2]);
% 
     [m,p_i]=min(abs(fr-f)) ;
     res(j,:)=[alpha_x(act_trap,p_i),alpha_y(act_trap,p_i)];
     res_all_x(j,:)=alpha_x(:,p_i)';
     res_all_y(j,:)=alpha_y(:,p_i)';
     f_out(j)=f;

end
loglog(f_out,abs(real(res_all_x(:,:))));
hold on
loglog(f_out,abs(imag(res_all_x(:,:))),'r');
hold off
