function plotBound(fig_switch,debug)

    if debug
        vSubject = {'A'};
        vStage = {'Pre'};
    else
        vSubject = {'A','B','C','N'};
        vStage = {'Pre','TimeLimit','Post'};
    end
    
    vBoundType = {'Optimal','Empirical',};
    vColor = {'k','r','b'};
    root = get_root;
    data_dir_nDDR = fullfile(root,'nDDM','results','empirical_max_likelihood','ndt_s_fit','params');
    
    switch fig_switch
        case 1
            saveBound;
        case 2
            plotBound;
    end

    function saveBound
        for si = 1:length(vSubject)
            for tti = 1:length(vStage) % 1:Pre, 2:TimeLimit, 3:Post
                clear Bup Blo
                % Plot optimal bounds
                data_dir_optim = fullfile(root,'OptimBound','results');
                file_name = sprintf('%s_%s.mat',vSubject{si},vStage{tti});
                a = dir(fullfile(data_dir_optim,file_name));
                fprintf('Loading %s created on %s\n',file_name,a.date);
                load(fullfile(data_dir_optim,file_name),'opt');
                DVDELTA_FINE = 0.02;
                MAXDV = 2;
                dv_axis_fine = -MAXDV:DVDELTA_FINE:MAXDV;            
                
                for ti = 1:size(opt.bestQ,2)
                    if any(opt.bestQ(:,ti)==3)
                        Blo(ti) = dv_axis_fine(find(opt.bestQ(:,ti)==2,1,'last'));
                        Bup(ti) = dv_axis_fine(find(opt.bestQ(:,ti)==1,1,'first'));
                    else
                        Blo(ti) = 0;
                        Bup(ti) = 0;
                    end
                end
                Bup = [nan,Bup];
                Blo = [nan,Blo];
                
                BUP.(vSubject{si}).(vStage{tti}).(vBoundType{1}) = Bup;
                BLO.(vSubject{si}).(vStage{tti}).(vBoundType{1}) = Blo;

                maxT = 10;
                dt = 0.0005;
                t_axis = 0:dt:maxT;
                
                b = dir(fullfile(data_dir_nDDR,file_name));
                fprintf('Loading %s created on %s\n',file_name,b.date);
                load(fullfile(data_dir_nDDR,file_name),'P')

                BUP.(vSubject{si}).(vStage{tti}).(vBoundType{2}) = P.Bup(1:length(t_axis));
                BLO.(vSubject{si}).(vStage{tti}).(vBoundType{2}) = P.Blo(1:length(t_axis));
                
                optim.(vSubject{si}).(vStage{tti}) = P.optim;
                
                fprintf('%s %-10s kappa = %.2f ndt_m = %.2f ndt_s = %.4f fval = %.4f\n',...
                    vSubject{si},vStage{tti},P.optim.theta(1),P.optim.theta(2),P.optim.theta(3),P.optim.fval);
            end
            fprintf('\n\n');
        end
        file_name = 'subject_bound.mat';
        save(fullfile(data_dir_nDDR,file_name),'BUP','BLO','t_axis','optim');

    end

    function plotBound
        % vSubject = {'A','B','C','N',};
        n_sub = 4;
        [p,~,~] = fileparts(mfilename('fullpath'));
        file_name = 'subject_bound.mat';
        load(fullfile(data_dir_nDDR,file_name),'BUP','BLO','t_axis','optim');
        
        dt = t_axis(2)-t_axis(1);
        
        figure(1);clf;

        for si = 1:length(vSubject)
            for ti = 1:length(vStage) % 1:Pre, 2:TimeLimit, 3:Post
                
                % find an appropriate time range for plotting (1-99th pecentile of RT)
                subject = vSubject{si};
                stage   = vStage{ti};
                pars.root = get_root;
                dat = make_data(subject,stage,[],pars);
                rt_lo = prctile(dat.rt,1);  %  1st percentile
                rt_hi = prctile(dat.rt,99); % 99th percentile
                
                for bi = 1:2 % 1:Optimal, 2:Empirical
                    subplot(2,4,n_sub*(bi-1)+si);hold on
                    ls = '-';
                                            
                    Bup = BUP.(vSubject{si}).(vStage{ti}).(vBoundType{bi});
                    Blo = BLO.(vSubject{si}).(vStage{ti}).(vBoundType{bi});
                    nt = length(Bup);
                    t_ax = (1:nt) .* dt;
                    % pick t_ax for an appropriate range
                    pick = false(size(t_ax));
                    pick(1:20:end) = true; % down-sample
                    pick = pick & t_ax>rt_lo & t_ax<rt_hi; % pick t_ax for an appropriate range
                    plot(t_ax(pick),Bup(pick),'LineStyle',ls,'color',vColor{ti},'LineWidth',1);
                    plot(t_ax(pick),Blo(pick),'LineStyle',ls,'color',vColor{ti},'LineWidth',1);
                    
                    ytick = -1.5:0.5:1.5;
                    labels = string(ytick); % extract
                    labels(1:2:end) = nan;
                    set(gca,'YTick',ytick,'YTickLabel',labels);
                    if si>1
                        set(gca,'YTickLabel',[]);
                    end
                    
                    % show logL in the upper right corner
                    txt = sprintf('fval = %.4f',optim.(vSubject{si}).(vStage{ti}).fval);
                    xpos = 2;
                    ypos = 1.5-0.25*ti;
                    
                    xlim([0 3]);
                    ylim([-1.5 1.5]);
                    set(gca,'Box',0)
                    
                    if si==1
                        switch bi
                            case 1
                                ylabel('Optimal')
                            case 2
                                ylabel('Empirical (nDDM)')
                        end
                    end
                    if bi==1
                        title(vSubject{si});
                    end
                end
                
                opt = optim.(vSubject{si}).(vStage{ti});
                
                fprintf('%s %-10s kappa = %.2f ndt_m = %.2f ndt_s = %.3f fval = %.4f\n',...
                    vSubject{si},vStage{ti},opt.theta(1),opt.theta(2),opt.theta(3),opt.fval);
                
            end
            fprintf('\n\n');
                        
        end
        set(gcf,'OuterPosition',[100 500 1050 450])
    end

end
