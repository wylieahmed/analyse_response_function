function [kappa,sigma]=load_response_data(ddir)


rsquare_thresh=0.998
fs=22;
eta_i=0.014;
eta_o=0.001;
s=1e-8;
k=3e-18;
r=3.5E-6;
cutoff=1000;
cutoff_low=.1;
fit_array=[1 1 1 0]
lm=round(r/.7e-7);
%Here I will oad all the  *deformation_response.mat files
% Which are supposed to contain the info required to reconstruct the
% Response with new calibration values.

 if (nargin == 0)
     %[ddir]=uigetdir('E:\Science\data\response_function');
    ddir='E:\Science\data\response_function\rbc_2µm\rbc_2µm_101007_lot_of_noise\1\response_function_series_1005'
 end
files=dir([ddir,filesep,'*deformation_response.mat']);

%now I open each of the files and store the important data in a structure
%callred res_data
for j=1:length(files)
    load([ddir,filesep,files(j).name])
    res_data(j).data=data;
    res_data(j).scan_path=scan_path;
    %res_data(j).repead_oss=repead_oss;
    res_data(j).f=f;
    res_data(j).a_x=a_x;
    res_data(j).a_y=a_y;
    res_data(j).s_eff=s_eff;
end


%new I will read the results tdms file
[rdata]=convertTDMS(0,[ddir,filesep,'results.tdms']);

%in(1,:)=rdata.Data.MeasuredData(15).Data;
%in(2,:)=rdata.Data.MeasuredData(16).Data;



for i=1:length(res_data)
    f(i)=res_data(i).f;
    a_x(i)=res_data(i).a_x(1);
    a_y(i)=res_data(i).a_y(1);
end

%sort the data
s_a(1,:)=f;
s_a(2,:)=a_x;
s_a(3,:)=a_y;

s_s=sortrows(s_a')';
f=s_s(1,:);
a_x=s_s(2,:);
a_y=s_s(3,:);


load('E:\Programing\Matlab\analyze_response_function\normal_rbc.mat')
t(:,2:3)=t(:,2:3);%*1e-18
[sigmat, kappat,eta, l_m, xdatat,ydata,FittedCurvet] = fit_psd_spherical_harmonics_safran_variable_parameters(t(:,1)',t(:,2)', r, s,k,.5*(eta_o+eta_i),lm,cutoff,cutoff_low,fit_array)


psd_res=4e-20./(pi*f).*abs(imag(a_x))*.01;
%now I tru to fit the PSD from the response data
[sigma, kappa,eta, l_m, xdata,ydata,FittedCurve] = fit_psd_spherical_harmonics_safran_variable_parameters(f,psd_res*1e18, r, s,k,.5*(eta_o+eta_i),lm,cutoff,cutoff_low,fit_array)
[f_out,FittedCurve]=generate_powerspectrum_spherical_harmonics_safran(r,xdatat,kappa,sigma,eta,l_m);

hold off
loglog(xdatat,FittedCurvet*1e18,'-k','LineWidth',2);
hold on
loglog(f_out,FittedCurve*1e18,'-r',f,psd_res*1e18,'*r','LineWidth',2);%,f,4e-20./(pi*f).*abs(imag(a_y)));
errorxy(t,'MarkSize',6,'EdgeEB',0.0,'Marker','d','EdgeColor','k','ScaleX','log','ScaleY','log')


xlabel('Frequency [Hz]');
ylabel('PSD [nm^2/Hz]');
legend('PSD direct measure','PSD from \alpha')
        set(gca,'FontSize',fs+5)
        set(gca,'LineWidth',2)
        set(gca,'XLim',[0.1 1000])
        set(gca,'YLim',[8E-4 3E3])
        set(gca,'XTick',[.1 1 10 100 1000])
        set(gca,'YTick',[.01 1 100])
        %title(['Normal RBC, \sigma=',num2str(abs(sigma)),' , \eta= ',num2str(abs(eta)),' , \kappa=',num2str(abs(kappa)),' , l_{max}=',num2str(r/l_m)])
        title('Hints for violation of the FDT in RBC')


%now I tru to fit the PSD from the response data

sigma
