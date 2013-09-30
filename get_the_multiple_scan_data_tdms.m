function [fm,pm]=get_the_multiple_scan_data_tdms(path)

path='E:\Science\data\response_function\3.28mu_beads_110617\6_same_as3_but_new_bead\multiple_run_series_1000'
%path='E:\Science\data\response_function\3.28mu_beads_110617\5_same_as3_but_after_lost_cell\multiple_run_series_1000'
%path='E:\Science\data\response_function\3.28mu_beads_110617\4_same_as3_but_after_lost_cell\multiple_run_series_1000'
%path='E:\Science\data\response_function\3.28mu_beads_110617\3\multiple_run_series_1011'
%path=uigetdir('E:\Science\data\response_function\3.28mu_beads_110617\','Please specify the parent dir');

files_tdms=dir([path,'\fluctuation_data*.tdms']);
files_lvb=dir([path,'\data_str*.lvb']);

pxm=0;
pym=0;
for i=1:length(files_tdms)
    [rdata]=convertTDMS(0,[path,'\',files_tdms(i).name]);
    eff_s= rdata.Data.Parameters.Root.Effective_Sampling_Rate.value
    [anum,param,data_r] = lbv_bin_read_with_parameters([path,'\',files_lvb(i).name]);
    l=length(rdata.Data.MeasuredData)/4
    for j=1:l
        in(j,:)=rdata.Data.MeasuredData(j).Data;
          if (in(j,3)==255)
             slo=rdata.Data.MeasuredData(4+3*j).Data
            [f_p,px]=power_sd(data_r(1,j:l:end)./data_r(3,j:l:end)/(slo(1)*1e6),eff_s/l);
            [f_p,py]=power_sd(data_r(2,j:l:end)./data_r(3,j:l:end)/(slo(2)*1e6),eff_s/l);
            [fl,pxl]=log_binning(f_p,px,200);
            [fl,pyl]=log_binning(f_p,py,200);
          end
    end
    pxm=pxm+pxl;
     pym=pym+pyl;
 
    hold on
    
end
   loglog(fl,pxl./length(files_tdms),fl,pyl./length(files_tdms))
hold off
