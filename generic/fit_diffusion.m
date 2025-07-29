function fit_diffusion(fn_wrapper,tl,th,tg,debug)

paralle_flag = 1; % use parallel toolbox?

if debug
    vSubject = {'A'};
    vStage = {'Pre'};
else
    vSubject = {'A','B','C','N'};
    vStage = {'Pre','TimeLimit','Post'};
end
IND = combvec(1:length(vSubject),1:length(vStage));

subjects = vSubject(IND(:,1));
stages   = vStage(IND(:,2));
n = size(IND,1);

pars = make_pars(debug);
pars.debug = debug;

g = InitGuessClass(tl,th,tg);
g.generate_tguess_random(pars.initguess_random_N,pars.initguess_random_seed,n);

%% run the optim
% Use parfor if necessary
cd(pars.temp_dir);
if 1
    if debug || ~paralle_flag
        % Run normal for loop
        for i=1:g.N
            dat = make_data(subjects{g.group(i)},stages{g.group(i)},[],pars);
            opt = do_fit(fn_wrapper,dat,tl,th,g.TG(i,:),pars);
            saveX(g.temp_filename{i},opt);
        end
    else
        % Use parallel toolbox
        parfor i=1:g.N
            dat = make_data(subjects{g.group(i)},stages{g.group(i)},[],pars);
            opt = do_fit(fn_wrapper,dat,tl,th,g.TG(i,:),pars);
            saveX(g.temp_filename{i},opt);
        end
    end

    %% gather files 
    g.gather_files();
    g.delete_temp_files();
else
    g.get_optim_filename(n);
end

%% eval best and save
if 1
    for i=1:n
        dat = make_data(subjects{i},stages{i},[],pars);
        plot_flag = true;
        opt = load(g.gather_filename{i});
        theta = opt.mejor.theta;
        [~,P] = feval(fn_wrapper,theta,dat.rt,dat.coh,dat.choice,dat.c,pars,plot_flag);

        P.optim = opt.mejor;%solo por compatibilidad
        % save fig and data
        
        save_dir = '../figures'; makedir(save_dir);
        figname = fullfile(save_dir,['fig_',subjects{i},'_',stages{i}]);
        saveas(gcf,figname,'fig')
        print(gcf,figname,'-depsc')
        close(gcf);
        
        save_dir = '../params'; makedir(save_dir);
        filename = fullfile(save_dir,[subjects{i},'_',stages{i}]);
        save(filename,'P','opt');
    end
end

%% Compute smooth RT distributions and plot psych & chrono curves
run_calc_smooth(vSubject,vStage,pars);
run_plot_smooth_perf(vSubject,vStage,pars);

end

function saveX(filename,opt)
    %trick to save inside par loop
    theta = opt.theta;
    fval = opt.fval;
    exitflag = opt.exitflag;
    output = opt.output;
    tl = opt.tl;
    th = opt.th;
    tg = opt.tg;
    save(filename,'theta','fval','exitflag','output','tl','th','tg');
end


