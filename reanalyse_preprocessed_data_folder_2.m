function reanalyse_preprocessed_data_folder_2(path)

files=dir([path,'\*.mat']);

for i=1:length(files)
    load([path,'\',files(i).name]);
    for j=1:length(response)
        p=response(j).path;
        %[res,f]=process_response_data_from_raw_data(p)
        [res,f,res_all_x,res_all_y,xy_k_extract]=process_response_data_plus_trap_calibration(p)
%         loglog(f,abs(real(res(:,1))));
%         hold on
%         loglog(f,imag(res(:,1)),'r');
%         hold off
        response(j).res_all_x=res_all_x;
        response(j).res_all_y=res_all_y;
        response(j).alphax=res(:,1);
        response(j).alphay=res(:,2);
        response(j).f=f;
        response(j).xy_k_extract=xy_k_extract
        pause(0.1)
        save([path,'\',files(i).name],'fluctuation','o_slope','override','response');
        
    end
    clear response
end
    