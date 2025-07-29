function demo

    % This code generate simulated data with known ground truth parameters:
    % kappa (drift term coefficient), the mean and s.d. of non-decision
    % time, and the bound shape.
    % Then it uses Dynamic Programming (d2b_dp) to computes the optimal 
    % bound shape (optimal policy) that maximizes the reward rate.
    
    % If you have not compile a mex function for Dynamic Programming,
    % run the following script.
    
    % compile_d2b_dp;
    
    debug = 1; % If set to 1, the dy
    
    rng(0);
    
    subject = 'demo';
    stage = 'Pre'; % 'Pre' or 'TimeLimit'
    
    %%%%%%% Ground truth parameters used for Monte-Carlo simulation %%%%%%%
    kappa = 15;
    ndt_m = 0.2; % (s)
    ndt_s = 0.01; % (s)
    
    maxT = 5; % (s)
    dt = 0.0005; % (s)
    t_axis = dt:dt:maxT; % (s)
    tau = 3; % (s)
    Bup = 1.2*exp(-t_axis./tau); % upper bound
    Blo = -Bup; % lower bound (symmetric)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cohs_pos = [0,2.^(5:9)]/1e3;
    cohs = [fliplr(-cohs_pos),cohs_pos];
    
    
        
    % For deadlines (sampled from a Rayleigh distribution)
    fshift = 0.6; % (s)
    sigma = 1.0; % (s)  
    
    % Experimental settings
    t_iti = 2; % (s) inter-trial interval
    t_misc_val = 0.7; % (s) miscellaneous time for stimulus latencies and fixation on valid (non-canceled) trials
    t_misc_inv = 0.4; % (s) miscellaneous time on invalid (canceled) trials
    t_iti_misc = t_iti + [t_misc_val, t_misc_inv]; % (s)    
    
    %%%%%%%%%%%%%%%% Create data by Monte-Carlo simulation %%%%%%%%%%%%%%%%

    num_trial = 1000; % per coherence

    RT = nan(num_trial,length(cohs)/2);
    RT_cor = cell(1,length(cohs)/2);
    RT_err = cell(1,length(cohs)/2);
    RT_inv = cell(1,length(cohs)/2);
    result = nan(size(RT));
    correct = nan(size(RT));

    t_dl = false(num_trial,length(t_axis));
    t_dl(:,end) = true;

    %%

    for ci = 1:length(cohs)

        fprintf('simulating coh = %.1f%%\n\n',cohs(ci)*100)

        % initialize the variables
        up_chosen = false(num_trial,1);
        lo_chosen = false(num_trial,1);
        chosen = false(num_trial,1);
        aborted = false(num_trial,1);
        ongoing = false(num_trial,length(t_axis));

        temp = kappa*cohs(ci)*dt + sqrt(dt)*randn(size(ongoing));
        nt = maxT/dt;
        dv = cumsum(temp,2);
        for ti = 1:nt

            aborted_now = ~aborted & ~chosen & t_dl(:,ti);
            aborted = aborted | aborted_now;

            up_chosen_now = ~aborted & ~chosen & dv(:,ti)>Bup(ti);
            lo_chosen_now = ~aborted & ~chosen & dv(:,ti)<Blo(ti);

            chosen = chosen | up_chosen_now | lo_chosen_now;
            up_chosen = up_chosen | up_chosen_now;
            lo_chosen = lo_chosen | lo_chosen_now;

            ongoing(:,ti) = ~aborted & ~chosen & (dv(:,ti)<Bup(ti) & dv(:,ti)>Blo(ti));

        end
        DT  = nansum(ongoing,2)*dt; % decision time
        nDT = ndt_m + ndt_s.*randn(sum(num_trial),1); % non-decision time
        RT(:,ci) = DT + nDT;

        result(up_chosen,ci) = 1;
        result(lo_chosen,ci) = 0;
        result(aborted,ci)   = nan;
        if cohs(ci)~=0
            correct(:,ci) = result(:,ci)==(cohs(ci)>0);
        else
            correct(:,ci) = result(:,ci)==(rand(num_trial,1)>0);
        end

        RT_cor{ci} = RT(result==1);
        RT_err{ci} = RT(result==0);
        RT_inv{ci} = RT(isnan(result));

        P_cor(ci) = sum(correct(:,ci)==1)/num_trial;
        P_err(ci) = sum(correct(:,ci)==0)/num_trial;
        P_inv(ci) = sum(isnan(correct(:,ci)))/num_trial;
        P_chosen(ci) = (sum(up_chosen) + sum(lo_chosen))/num_trial;
        
    end

    RT_mean = nanmean(RT(:));
    Pcor_mean = mean(P_cor(:));
    Perr_mean = mean(P_err(:));
    PR_noDL = (Pcor_mean - Perr_mean)/(RT_mean + t_iti_misc(1));
    fprintf('Monte-Carlo simulation w/ deadline:\nPR = %f\n',PR_noDL)

    % Assemble data across coherences
    choice = []; coh = []; rt = []; c = [];
    for ci = 1:length(cohs)
        inds = ismember(result(:,ci),[0,1]) & ~isnan(RT(:,ci));
        choice = [choice; result(inds,ci)];
        coh    = [coh; cohs(ci).*ones(sum(inds),1)];
        rt     = [rt; RT(inds,ci)];
        c      = [c; correct(inds,ci)];
    end

    dat.choice = logical(choice);
    dat.coh    = coh;
    dat.rt     = rt;
    dat.c      = c;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % non-parametric DDM (nDDM) has three fitted parameters in theta.
    % theta(1): kappa
    % theta(2): mean of non-decision time
    % theta(3): s.d. of non-decision time

    % tl: lowest value for theata for fit
    % th: highest value for theta for fit

    % During the fit, tl and th provides the lower and upper bounds.
    % If they are the same value, the parameter is fixed at that value during 
    % the fit. See fminsearchbnd.m for details.
 
    tl = [  1, 0.001, 0.001, 0, 0];
    th = [100,     1,     1, 0, 0];
    tg = [];

    pars = make_pars(0);
    
    init_guess = [10, 0.1, 0.1, 0, 0]; % initial guess for the fit
    plot_flag = 0;
    
    % Fitting a nDDM.
    % Once you save the result, you can load the result from next time.
    if 0
        opt = do_fit(@wrapper_dtb_empiricalbound_rt, dat, tl, th, init_guess, pars);
        [~,P] = wrapper_dtb_empiricalbound_rt(fn_wrapper, opt.theta, dat.rt, dat.coh, dat.choice, dat.c, pars, plot_flag);
        P.optim = opt;
        save('opt.mat','P');
    else
        load('opt.mat','P');
        % [err,~] = feval(fn_wrapper,[kappa, ndt_m, ndt_s, 0, 0], dat.rt, dat.coh, dat.choice, dat.c, pars, plot_flag);
    end    
    rt_lo = prctile(dat.rt,1);  %  1st percentile
    rt_hi = prctile(dat.rt,99); % 99th percentile
    
    figure(2);clf;hold on;
    pick1 = t_axis>rt_lo & t_axis<rt_hi;
    pick2 = t_axis>rt_lo & t_axis<rt_hi;
    
    plot(t_axis(pick1),Bup(pick1),'b-')
    plot(P.t(pick2),P.Bup(pick2),'r-')
    plot(t_axis(pick1),Blo(pick1),'b-')    
    plot(P.t(pick2),P.Blo(pick2),'r-')
    axis([0 3 -1.5 1.5]);
    legend('Ground truth','fit by nDDM');
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Dynamic Programming %%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Running DP\n\n ')
    
    % Get fitted parameters
    kappa = P.optim.theta(1);
    ndt_m = P.optim.theta(2);
                
    if debug
        % coarse but faster
        MAXDV = 2;
        DVDELTA = 0.05;
        maxT = 5; % (s)
        dt = 0.01; % (s)
    else
        MAXDV = 2;
        DVDELTA = 0.02;
        maxT = 5; % (s)
        dt = 0.0005; % (s)
    end
    
    t_axis = dt:dt:maxT;
    dv_axis = -MAXDV:DVDELTA:MAXDV;
    coherences = cohs;
    p_dl = zeros(1,length(t_axis));

    rewardSuccess     = 1; % subjects get one point for a correct choice
    rewardFailure     = -1; % subjects lose one point for an error choice
    rewardListen      = 0; % cost of accumulation
    rewardTimesUp     = 0;
    gamma             = 1.0; % temporal discount factor (no discount when gamma = 1.0)
    maxError          = 0.00001;
    rewardPrm         = [rewardSuccess,rewardFailure,rewardListen,rewardTimesUp];
    
    t_iti_nd_misc_val = t_iti + ndt_m + [t_misc_val, t_misc_inv]; % Include ndt_m in ITI
    
    switch stage
        case 'Pre'
            pdf = zeros(size(t_axis));
            hr = zeros(size(t_axis));
            sf = ones(size(t_axis));
        case 'TimeLimit'
            % Dealdline is set on a half of trials.
            % Note that ndt_m is subtracted from tshift (t0) of the 
            % Rayleigh distribution.
            pdf = rayleighPDF(t_axis,dt,fshift - ndt_m, sigma)./2;
            sf = rayleighSF(t_axis,fshift - ndt_m, sigma)./2 + 0.5;
            hr = pdf./sf; % hazard rate
            % hr(end) = 1; % The trial terminates at 5 sec.
    end
        
    %%%%%%%%%%%%%%%%% Configure file name for TlistenVal %%%%%%%%%%%%%%%%%%
    file_name = sprintf('TlistenVal_%s_%s_MAXDV=%.2f_DVDELTA=%.2f_maxT=%.2f_dt=%.2f_cohs=',subject,stage,MAXDV,DVDELTA,maxT,dt);
    cohs_pos_int = cohs_pos*1e3;
    for ci = 1:length(cohs_pos_int)
        if ci==1
            file_name = sprintf('%s%d',file_name,cohs_pos_int(ci));
        else
            file_name = sprintf('%s_%d',file_name,cohs_pos_int(ci));
        end
    end
    file_name = [file_name,'.mat'];
    save_dir = fullfile(get_root,'nDDM\results\empirical_max_likelihood\ndt_s_fit','TlistenVal');
    makedir(save_dir);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if debug
        fn_d2b = @d2b_dp;
    else
        fn_d2b = @d2b_dp_mex; % run faster with MEX code 
    end
    
    % Creating a trnsition matrix (TlistenVal) takes long time.
    % If it does not exist, create and save it, and load the saved matrix from the next time.
    if exist(fullfile(save_dir,file_name),'file')
        fprintf('TlistenVal file found. Loading... ');
        load(fullfile(save_dir,file_name),'TlistenVal');
        fprintf('Done\n\n ')
    else
        TlistenVal = [];
        bestQ = [];
        temp = zeros(size(t_axis)); % place holder for hazard rate
        [~,~,~,~,~,TlistenVal] = feval(fn_d2b, MAXDV, DVDELTA, maxT, dt, coherences, kappa, gamma, rewardPrm, maxError, temp, t_iti_nd_misc_val, bestQ, 0);
        
        info.subject = subject;
        info.stage = stage;
        info.kappa = kappa;
        info.MAXDV = MAXDV;
        info.DVDELTA = DVDELTA;
        info.maxT = maxT;
        info.dt = dt;
        
        save(fullfile(save_dir,file_name),'TlistenVal','info','-v7.3');
    end
    
    TlistenVal = reshape(TlistenVal,size(TlistenVal,1),[]);
    bestQ = []; % compute the optimal policy
    
    t_iti_nd_misc_val = [t_iti + t_misc_val, t_iti + t_misc_inv] + ndt_m; % Include ndt_m in ITI
    
    [opt.bestQ,opt.Q,opt.V,opt.Vdelta,opt.rho,~] = feval(fn_d2b, MAXDV, DVDELTA, maxT, dt, coherences, kappa, gamma, rewardPrm, maxError, hr, t_iti_nd_misc_val, bestQ, TlistenVal);
    
    fprintf('Dynamic Programming predicts PR = %f\n',opt.rho)
    save_dir = fullfile(get_root,'OptimBound','results'); makedir(save_dir);
    file_name = sprintf('%s_%s.mat',subject,stage);
    save(fullfile(save_dir,file_name),'opt');
    
    if 1
        % plot the optimal policy
        figure(3);clf
        imagesc(t_axis,dv_axis,opt.bestQ,[1 3])
        axis xy
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
