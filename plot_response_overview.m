function plot_response_overview


path{1}='E:\Science\data\response_function\pre_analysed_data\healthy_cells';
path{2}='E:\Science\data\response_function\pre_analysed_data\starved_and_depleted_old';
path{3}='E:\Science\data\response_function\pre_analysed_data\starved_new';

for k=1:3
files=dir([path{k},'\*.mat']);
n=0;

    for i=1:length(files)
        load([path{k},'\',files(i).name]);
        
        for j=1:length(response)
            load([response(j).path,'\parameters.mat']);
            fi=dir([response(j).path,'\*deformation*.mat'])
            load([response(j).path,'\',fi(1).name])
            act=response(j).act_trap;
            n=n+1;
            col_f(n,:)=response(j).f;
            col_res(n,:)=response(j).alphax;
            col_res_corr(n,:)=response(j).alphax*xy_k(act,1)/response(j).xy_k_extract(act,1);
            %[res,f]=process_response_data_from_raw_data(p)
            loglog(col_f(n,:),imag(col_res(n,:)))
            hold on
            loglog(col_f(n,:),imag(col_res_corr(n,:)),'r')

        end
        clear response
    end
    f_out(k,:)=mean(col_f,1);
    alphax(k,:)=mean(col_res,1);
    alphax_std_im(k,:)=std(imag(col_res),1);
    alphax_std_real(k,:)=std(real(col_res),1);
    alphax_corr(k,:)=mean(col_res_corr,1);
    alphax_corr_std_im(k,:)=std(imag(col_res_corr),1);
    alphax_corr_std_real(k,:)=std(real(col_res_corr),1);
    
    clear col_f,col_res,col_res_corr;
end
hold off
loglog(f_out,imag(alphax));
hold on
loglog(f_out,imag(alphax_corr));
save('E:\Science\data\response_function\pre_analysed_data\collection_with_corrected_trap_stiffness.mat','f_out','alphax','alphax_corr','alphax_std_im','alphax_std_real','alphax_corr_std_im','alphax_corr_std_real')