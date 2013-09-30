function reanalyse_preprocessed_data_folder_3(path)
%here we reprocess the response and the fluctuation data
%path='E:\Science\data\response_function\pre_analysed_data\healthy_cells'
path='E:\Science\data\response_function\pre_analysed_data\starved_new'
%path='E:\Science\data\response_function\pre_analysed_data\starved_and_depleted_old'
files=dir([path,'\*.mat']);

ind=0

for i=1:length(files)
    load([path,'\',files(i).name]);

    for j=1:length(response)
        p=response(j).path;
        %[res,f]=process_response_data_from_raw_data(p)
        f_fit=10;
        [res,f,res_all_x,res_all_y,xy_k_extract,fbl,Pbl,k_fit]=process_response_data_plus_trap_calibration(p,f_fit)
%         loglog(f,abs(real(res(:,1))));
%         hold on
%         loglog(f,imag(res(:,1)),'r');
%         hold off
        response(j).k_fit=k_fit
        response(j).res_all_x=res_all_x;
        response(j).res_all_y=res_all_y;
        response(j).alphax=res(:,1);
        response(j).alphay=res(:,2);
        response(j).f=f;
        response(j).xy_k_extract=xy_k_extract
        response(j).fbl=fbl;
        response(j).Pbl=Pbl;
        pause(0.1)
        save([path,'\',files(i).name],'fluctuation','o_slope','override','response');
        
    end
    
    for j=1:length(fluctuation)
        %allignment_frequency
        
        fa=40
        if strcmp(path,'E:\Science\data\response_function\pre_analysed_data\starved_new')
            fa=120
        end
        
        p=fluctuation(j).path;
        d=fluctuation(j).data;
        for k=1:length(response)
            [m_val,i_min]=min(abs(fa-response(k).fbl(response(k).act_trap,:)));
            c_corr_i(k)=response(k).Pbl(response(k).act_trap,i_min)*1e-12;
            cf_corr_i(k)=response(k).fbl(response(k).act_trap,i_min);
        end
        at=response(k).act_trap
        c_corr=mean(c_corr_i);
        cf_corr_i=mean(cf_corr_i);
        ind=0
        for k=1:length(d)
            ind=ind+1;
            [m_val,i_min]=min(abs(fa-d(k).f));
            f_out(ind,:)=d(k).f;
            px(ind,:)=d(k).px;
            slope(ind)=d(k).slopes(1);
            px_corr(ind,:)=d(k).px*(c_corr/px(ind,i_min));
        end
        fm=mean(f_out);
        pxm=mean(px);
        clear f_out
        clear px
       
        
        
        hold off
        loglog(fluctuation(j).f,fluctuation(j).psdxy(1,:),response(1).fbl(at,:)',response(1).Pbl(at,:)'*1e-12)
        hold on
        for i2=1:length(response)
            loglog(response(i2).f,imag(response(i2).alphax)*4.1e-21./(response(i2).f'*pi),'r')
        end        
        [m_val,i_min]=min(abs(fa-fluctuation(j).f))
        fluctuation(j).psd_corrected=fluctuation(j).psdxy(1,:)*(c_corr/fluctuation(j).psdxy(1,i_min))
        
        %now I combine the high frequency from the high laser power with
        %the low frequency from the low laser power. The overlap  frequency
        %region will be 50-100 Hz
        f_min=40;
        f_max=80;
        if strcmp(path,'E:\Science\data\response_function\pre_analysed_data\starved_new')
            f_min=60;
            f_max=200;
        end
        [m_val,i_min]=min(abs(f_max-fluctuation(j).f))
        fluctuation(j).f_low=fluctuation(j).f(1:i_min);
        fluctuation(j).psd_low=fluctuation(j).psd_corrected(1:i_min);
        [m_val,i_min]=min(abs(f_min-response(1).fbl(at,:)));
        fluctuation(j).f_high=response(1).fbl(at,i_min:end)';
        fluctuation(j).psd_high=response(1).Pbl(at,i_min:end)'*1e-12;
        hold off
        
        %here I create a sorted psd anf frequency vector that will be used
        %for the paper. This then extends the data to the high frequencies
        %and it is justified by the fact that the high frequency data for
        %strong and weak traps will overlap. This overlap was also used to
        %adjust the low laser data for possible noise effects on the
        %calibration
        psd_int(1,:)=[fluctuation(j).f_low fluctuation(j).f_high'];
        psd_int(2,:)=[fluctuation(j).psd_low fluctuation(j).psd_high'];
        psd_out=(sortrows(psd_int',1))'
        fluctuation(j).psd_final=psd_out;
        
        loglog(fluctuation(j).f_low,fluctuation(j).psd_low,'+r')
        hold on
        loglog(fluctuation(j).f_high,fluctuation(j).psd_high,'+k')
        pause(1)
    end
    save([path,'\',files(i).name],'fluctuation','o_slope','override','response');
    clear response
    clear psd_int;
end
    