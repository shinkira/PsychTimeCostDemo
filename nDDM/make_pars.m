function pars = make_pars(debug)

    pars.dt_flag = 1; %
    pars.error_flag = 2; %1 is perf, 2 is likelihood

    % pars.initguess_flag = 1;%1: use tg; 2: use grid; 3 use random
    pars.initguess_random_seed = 123456; %124321;
    pars.initguess_random_N = 10;
    model_dir = fullfile(get_root,'nDDM/results/empirical_max_likelihood/ndt_s_fit');
    temp_dir = fullfile(model_dir,'temp');
    makedir(temp_dir);
    pars.model_dir = model_dir;
    pars.temp_dir = temp_dir;
    pars.root = get_root;
    pars.debug = debug;
    
    if debug
        pars.initguess_random_N = 1;
    end

end