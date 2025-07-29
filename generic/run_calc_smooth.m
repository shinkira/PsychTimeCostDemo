function run_calc_smooth(vSubject,vStage,pars)

% vSubject = {'A','B','C','N'};
% vStage = {'Pre','TimeLimit','Post'};

for i=1:length(vSubject)
    subject = vSubject{i};
    for j=1:length(vStage)
        
        stage = vStage{j};
        
        %         dat = make_data(subject,stage);
        %         rt = dat.rt;
        %         coh = dat.coh;
        %         choice = dat.choice;
        
        file_name = sprintf('%s_%s.mat',subject,stage);
        file_path = fullfile(pars.model_dir,'params',file_name);
        load(file_path)
        
        %%
        kappa = P.optim.theta(1);
        ndt_m = P.optim.theta(2);
        coh0 = P.optim.theta(4);
        out = calc_smooth_perf_rt(P,kappa,ndt_m,coh0,pars);
        ucoh = out.ucoh;
        pright = out.pright;
        rt = out.rt;
        rt_c = out.rt_c;
        rt_nc = out.rt_nc;
        
        %% save fig and data
        save_dir = fullfile(pars.model_dir,'smooth'); makedir(save_dir);
        file_path = fullfile(save_dir,file_name);
        save(file_path,'ucoh','pright','rt','rt_c','rt_nc');
        
    end
end




end



