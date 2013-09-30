function [res,f_out,res_all_x,res_all_y,xy_k_extr,fbl,Pbl,K_fit]=process_response_data_plus_trap_calibration(path,f_fit)

%This function will open the raw data of the response function experiments,
%find the active trap, and try to get the opposing bead. Then get the real and the imaginary part of the active part and of the correlation between the other beads. 
%It will take the trap stiffness and the position calibration as found in
%the path. Output is the response structure as will be used later in
%analysis.
%Furthermore it will output a response structure that gives the correlation
%mixes between the different beads with the active bead. The datasets are
%either the absolute value of the active bead, or the distance values with
%respect to that bead
%This program will also take the fluctuation data and collect it to have a
%final psd of the free fluctuations that can be used to determine the trap
%stiffness. 

%path='E:\Science\data\response_function\3.28um_beads_rbc\3.28mu_beads_110726\1\response_function_series_1000'
files=dir([path,'\*response.mat']);
load([path,'\parameters.mat']);
act_trap=find(Parameter.Traps(:,3)==255)
n_traps=size(Parameter.Traps,1);
AODmu=0.0265;    
psd_x_coll=[];
psd_y_coll=[];
f_coll=[];


for j=1:length(files)
    %load the data
    load([path,'\',files(j).name]);
    %here I extract the trap and the bead positions of each of the traps
    clear x_trap y_trap power_trap bead_rel_x bead_rel_y abs_bead_x abs_bead_y Fx Fy alpha_x alpha_y 
        psd_x=[];
        psd_y=[];
    for i=1:n_traps

        x_trap(i,:)=squeeze(data(i,1,:))/AODmu;
        y_trap(i,:)=squeeze(data(i,2,:))/AODmu;
        power_trap(i,:)=squeeze(data(i,3,:));
        bead_rel_x(i,:)=-1*squeeze(data(i,4,:))./squeeze(data(i,8,:))/xy_slope(i,1);
        bead_rel_y(i,:)=-1*squeeze(data(i,6,:))./squeeze(data(i,8,:))/xy_slope(i,2);
        abs_bead_x(i,:)=x_trap(i,:)+bead_rel_x(i,:);
        abs_bead_y(i,:)=y_trap(i,:)+bead_rel_y(i,:);
        Fx(i,:)=1*xy_k(act_trap,1)*bead_rel_x(i,:);
        Fy(i,:)=1*xy_k(act_trap,2)*bead_rel_y(i,:); 
        
        [fp,psd_x(i,:)]=power_sd(bead_rel_x(i,:),s_eff);
        [fp,psd_y(i,:)]=power_sd(bead_rel_y(i,:),s_eff);
        %[m,p_i]=min(abs(fp-f)) ;

        
    end
    [mp,p_i]=min(abs(fp-f));
    psd_x(:,p_i)=[];
    psd_y(:,p_i)=[];
    f_int=fp;
    f_int(p_i)=[];
    if (f>f_fit & f<50)
        psd_x_coll=[psd_x_coll,psd_x(:,2:end)];
        psd_y_coll=[psd_y_coll,psd_y(:,2:end)];
        f_coll=[f_coll,f_int(2:end)];
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
     
     loglog(f_int,psd_x(3,:),'+')

end
pox(1,:)=f_coll;
pox(2:5,:)=psd_x_coll(1:4,:);
po_sx=sortrows(pox',1);
poy(1,:)=f_coll;
poy(2:5,:)=psd_y_coll(1:4,:);
po_sy=sortrows(poy',1);
sta=3;
for i=1:4
    [fb,Pb]=lin_binning(po_sx(:,1),po_sx(:,i+1),10000);
    xy_k_extr(i,1)=4e-21/(sum(diff(fb(sta:end)).*Pb(sta:end-1))*1e-12*2);
    [fbl(i,:),Pbl(i,:)]=log_binning(po_sx(:,1),po_sx(:,i+1),100);
    [fb,Pb]=lin_binning(po_sy(:,1),po_sy(:,i+1),10000);
    xy_k_extr(i,2)=4e-21/(sum(diff(fb(sta:end)).*Pb(sta:end-1))*1e-12*2);
end



%fit
[fb,Pb]=lin_binning(po_sx(:,1),po_sx(:,act_trap+1),3000);
logp=log(Pb);
options = fitoptions('Method','NonlinearLeastSquares');                     %fit options 
options.MaxIter = 10000;                                                    %options
options.StartPoint = [4.24e-3 10];                                        %options
model = fittype('log(D/(2*pi^2*(f^2+fc^2)))','indep','f');                    %model to fit
[c,gof] = fit(fb(1:100)',logp(1:100)',model,options);  %fit outputs

loglog(fb,Pb,'r+',fb,exp(c(fb)))

K_fit = 2*pi*3*pi*0.0009*3.28e-6*c.fc;
    


