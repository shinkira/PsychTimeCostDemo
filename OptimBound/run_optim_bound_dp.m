function run_optim_bound_dp(debug)

    % First, compute the point rate for a flat bound case with analytical
    % solutions. Then use the DP code to see what point rate it returns.
    
    % Deadline is set on a half of trials.
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
    
    if debug
        for i = 1:n
            runEachSubjectAndStage(subjects{i},stages{i},debug);
        end
    else        
        parfor i = 1:n
            runEachSubjectAndStage(subjects{i},stages{i},debug);
        end
    end

end

function runEachSubjectAndStage(subject,stage,debug)
        
    % get kappa
    [p,~,~] = fileparts(mfilename('fullpath'));
    data_dir = fullfile(p,'..','nDDM','results','empirical_max_likelihood','ndt_s_fit','params');
    file_name_1 = sprintf('%s_%s.mat',subject,stage);
    load(fullfile(data_dir,file_name_1),'P')
    
    kappa = P.optim.theta(1);
    ndt_m = P.optim.theta(2);
    ndt_s = P.optim.theta(3);
        
    % get ITI (including time from FP onset to the end of Tnd_mean)
    data_dir = fullfile(p,'..','SubjectData');
    file_name_2 = sprintf('%sData%s.mat',subject,stage);
    load(fullfile(data_dir,file_name_2),'data');
    [non_rt_iti_val, non_rt_iti_inv, ~] = calc_non_rt_iti(data);
    
    if 1
        fprintf('%s %s: ITI_val = %.2f ITI_inv = %.2f diff = %.3f\n',subject, stage, non_rt_iti_val, non_rt_iti_inv, non_rt_iti_val-non_rt_iti_inv);
        % return
    end
    
    if isnan(non_rt_iti_inv)
        % Avoid uisng NaN as ITI
        non_rt_iti_inv = non_rt_iti_val;
    end
    ITI = [non_rt_iti_val, non_rt_iti_inv] + ndt_m; % Include ndt_m in ITI 22/01/31

    cohs_pos = [0,2.^(5:9)]/1e3;
    cohs = [fliplr(-cohs_pos),cohs_pos];
    coherences = cohs;
    uv = kappa*cohs;
    
    pt_rate = get_observed_pt_rate(data);
    fprintf('pt_rate = %.4f %s %.4f\n',mean(pt_rate),char(177),std(pt_rate))
    
    % DP parameters
    MAXDV = 2;
    maxT = 5;
        
    if debug
        % coarse but faster (for debug)
        DVDELTA = 0.1;
        dt = 0.1;
    else
        DVDELTA = 0.02;
        dt = 0.0005; % 0.0001
    end
    
    t_axis = dt:dt:maxT;
    dv_axis = -MAXDV:DVDELTA:MAXDV;
    coherences = cohs;
    hr = zeros(1,length(t_axis));

    rewardSuccess     = 1;
    rewardFailure     = -1;
    rewardListen      = 0; % cost of accumulation
    rewardTimesUp     = 0;
    gamma             = 1.0;
    maxError          = 0.00001;
    rewardPrm         = [rewardSuccess,rewardFailure,rewardListen,rewardTimesUp];

    switch subject
        case 'C'
            fshift = 0.6; 
            sigma = 1.0;
        case 'N'
            fshift = 0.6;
            sigma = 0.6;
        case 'B'
            fshift = 0.3;
            sigma = 0.6;
        case 'A'
            fshift = 0.6;
            sigma = 0.7;
    end
    % ndt_m elapses before the diffusion process starts
    % thus subtract ndt_m to account for this.
    switch stage
        case 'Pre'
            pdf = zeros(size(t_axis));
            hr = zeros(size(t_axis));
            sf = ones(size(t_axis));
            
        case 'TimeLimit'
            if 0
                pdf = rayleighPDF(t_axis,dt,fshift - ndt_m, sigma);
                sf = rayleighSF(t_axis,fshift - ndt_m, sigma);
                hr = pdf./sf; % hazard rate
            else
                % Dealdline is set on a half of trials.
                % Note that ndt_m is subtracted from tshift (t0) of the 
                % Rayleigh distribution.
                pdf = rayleighPDF(t_axis,dt,fshift - ndt_m, sigma)./2;
                sf = rayleighSF(t_axis,fshift - ndt_m, sigma)./2 + 0.5;
                hr = pdf./sf; % hazard rate
                % hr(end) = 1; % The trial terminates at 5 sec.
            end
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
    get_root;
    save_dir = fullfile(get_root,'nDDM\results\empirical_max_likelihood\ndt_s_fit','TlistenVal');
    makedir(save_dir);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if exist(fullfile(save_dir,file_name),'file')
        fprintf('TlistenVal file found. Loading... ');
        load(fullfile(save_dir,file_name),'TlistenVal');
        fprintf('Done\n\n ')
        TlistenVal = reshape(TlistenVal,size(TlistenVal,1),[]);
    else
        TlistenVal = [];
        bestQ = [];
        temp = zeros(size(t_axis)); % place holder for hazard rate
        [~,~,~,~,~,TlistenVal] = d2b_dp_mex(MAXDV,DVDELTA,maxT,dt,coherences,kappa,gamma,rewardPrm,maxError,temp,ITI,bestQ,0);
        
        info.subject = subject;
        info.stage = stage;
        info.kappa = kappa;
        info.MAXDV = MAXDV;
        info.DVDELTA = DVDELTA;
        info.maxT = maxT;
        info.dt = dt;
        
        save(fullfile(save_dir,file_name),'TlistenVal','info','-v7.3');
    end
    
    fprintf('Running DP\n\n ')
    bestQ = []; % compute the optimal policy
    % Include ndt_m in ITI 22/01/31
    [opt.bestQ,opt.Q,opt.V,opt.Vdelta,opt.rho,~] = d2b_dp(MAXDV,DVDELTA,maxT,dt,coherences,kappa,gamma,rewardPrm,maxError,hr,ITI,bestQ,TlistenVal);
    
    fprintf('Dynamic Programming predicts PR = %f\n',opt.rho)
    save_dir = fullfile(get_root,'OptimBound','results'); makedir(save_dir);
    save(fullfile(save_dir,file_name_1),'opt');
    
    if 1
        % plot and save a fig of the optimal policy
        figure(1);clf
        imagesc(t_axis,dv_axis,opt.bestQ,[1 3]);colorbar
        axis xy
        fig_path = fullfile(get_root,'OptimBound','figures'); makedir(fig_path);
        fig_name = sprintf('%s_%s.mat',subject,stage);
        saveFig('png',fig_path,fig_name);        
    end

end

function pt_rate = get_observed_pt_rate(data)
    % compute experimental point rate
    pt = data(:,7);
    t_start = data(:,5);
    t_end = data(:,6);
    ind_start = [1;find(abs(diff(pt))>1)+1]; % index for the first trial
    ind_end   = [find(abs(diff(pt))>1);length(pt)];   % index for the last trial
    n_session = length(ind_start);
    for si = 1:n_session
        pt_earned = pt(ind_end(si)) - 50;
        t_session = t_end(ind_end(si)) - t_start(ind_start(si));
        pt_rate(si) = pt_earned/t_session;
    end
    % data
    % 1: result (1:correct, 0 error)
    % 2: choice direction (1:right, 2:left)
    % 3: coherence (pos: right, neg:left)
    % 4: RT (ms)
    % 5: time of trial start
    % 6: time of trial end
    % 7: point (start from 50)
    % 8: deadline (s)
end

function saveFig(format,fig_path,fig_name)
    
    h = findobj('type','figure');
    n = length(h);
    fprintf('Saving %s in %s...\n',fig_name,format);
    for i = 1:n
        figure(h(i).Number)
        h(i).Renderer='Painters';
        set(gcf,'color','w')
        switch format
            case 'png'
                fig_name = [fig_name,'.png'];
                print(fullfile(fig_path,fig_name),'-dpng');
            case 'pdf'
                fig_name = [fig_name,'.pdf'];
                print(fullfile(fig_path,fig_name),'-dpdf');
        end
    end

end