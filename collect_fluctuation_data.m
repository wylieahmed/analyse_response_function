path='E:\Science\data\RBC_blebbistatin\50_mue_mol\110909';
path='E:\Science\data\response_function\Kalpana\2013-09-17\bead2\';
dirs2=dir([path,filesep,'multiple_run_series*'])

    for j=1:length(dirs2)
        n=n+1;
        [f,psdxy,psdxy_std,slopes,act_trap]=process_fluctuation_folder_v3([path,'\',dirs2(j).name]);
        px(n,:)=psdxy(1,:);
        py(n,:)=psdxy(2,:);
    end

p_mean=mean(px,1)
loglog(f,p_mean,'k')