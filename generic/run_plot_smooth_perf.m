function p = run_plot_smooth_perf(vSubject,vStage,pars)

% vSubject = {'A','B','C','N'};
% vStage = {'Pre','TimeLimit'};

if exist('subjects','var')
    sub_ind = [];
    for i = 1:length(subjects)
        sub_ind = [sub_ind, find(strcmp(vSubject,subjects(i)))];
    end
else
    sub_ind = 1:length(vSubject);
end

if exist('stages','var')
    stg_ind = [];
    for i = 1:length(stages)
        stg_ind = [stg_ind, find(strcmp(vStage,stages(i)))];
    end
else
    stg_ind = 1:length(vStage);
end

for mi = 1

    switch mi
        case 1 % empirical bound
            figure(1);clf;
            color_set = [0 0 0;0 0 0];
            lw = 0.5;
            
            color = 'krb';
            color_nc = 0.4*[0,0,0;1,0,0;0 0 1] + 0.6*ones(3,3);
            
        case 2 % flat bound
            color_set = [0 0 1;0 0 1];
            lw = 1;
            
            color = 'ggg';
            color_nc = color_set;
    end    
    
    fh = figure(1);
    p = publish_plot(2,4,'hfig',fh);

    min_coh = 0.01;
    zero_coh = 0.005;

    for i = sub_ind
        
        subject = vSubject{i};
        
        for j = stg_ind

            stage = vStage{j};

            dat = make_data(subject,stage,[],pars);
            RT{i,j} = dat.rt;
            rt = dat.rt;
            coh = dat.coh;
            choice = dat.choice;
            c = dat.c;
            file_name = sprintf('%s_%s.mat',subject,stage);
            file_path = fullfile(get_root,'MiniData',file_name);
            load(file_path)

            % Load predicted mean RT and p_correct
            file_path = fullfile(pars.model_dir,'smooth',file_name);
            aux = load(file_path);
            ucohs = aux.ucoh;
            prights = aux.pright;
            rts = aux.rt;
            rts_c = aux.rt_c;
            rts_nc = aux.rt_nc;
            
            % Plot a psychometric curve
            p.current_ax(1,i);hold on
            tt = ucohs(ucohs>0);
            xx = 0.5*(prights(ucohs>0)+prights(ucohs(end:-1:1)<0));
            plot(tt(tt>min_coh),xx(tt>min_coh),'color',color(j),'LineWidth',lw);
            [tt,xx,ss] = curva_media(c==1,abs(coh),[],0); % plot predicted p_correct
            if mi==1
                % plot observed p_correct
                h = errorbar_noh(tt,xx,ss);
                set(h,'color',color(j))
                set(h(2),'LineStyle','none','Marker','.','MarkerSize',25,'color',color(j))
            end
            set(gca,'xscale','log')
            title(vSubject{i})
            
            % Plot a chronometric curve for correct trials
            p.current_ax(2,i);hold on
            xx = rts(ucohs==0);
            if mi==1
                plot(zero_coh,xx,'s','color',color(j)); % plot predicted RT at 0% coh
            end
            tt = ucohs(ucohs>0);
            xx = 0.5*(rts_c(ucohs>0)+rts_c(ucohs(end:-1:1)<0));
            plot(tt(tt>min_coh),xx(tt>min_coh),'color',color(j),'LineWidth',lw); % plot predicted RT for >0% coh
            [tt,xx,ss] = curva_media(rt,abs(coh),c==1,0);
            XX{i,j} = xx; % store mean RT across coherences
            rt0(i,j,1) = xx(1); % store 0% coh RT for y-axis adjustment
            tt(tt==0) = zero_coh;
            if mi==1
                h = errorbar_noh(tt,xx,ss); % plot observed RT
                set(h,'color',color(j))
                set(h(2),'LineStyle','none','Marker','.','MarkerSize',25,'color',color(j))
            end
            
            % Plot a chronometric curve for error trials
            unique_cohs = [0,2.^(5:9),inf]/1e3;
            hc = histcounts(abs(coh(c==0)),[0,2.^(5:9),inf]/1e3);
            coh_include = unique_cohs(hc>=5); % do not plot for cohs with less than n trials
            [tt,xx,ss] = curva_media(rt,abs(coh),c==0 & ismember(abs(coh),coh_include),0);
            rt0(i,j,2) = xx(1); % store 0% coh RT for y-axis adjustment
            tt(tt==0) = zero_coh;
            if mi==1
                h = errorbar_noh(tt,xx,ss);
            end

            tt = ucohs(ucohs>0);
            xx = 0.5*(rts_nc(ucohs>0)+rts_nc(ucohs(end:-1:1)<0));
            pick = tt>min_coh & tt<=max(coh_include)*(2^0.2);% Stop extrapolate too much beyond the data points
            plot(tt(pick),xx(pick),'color',color_nc(j,:),'LineWidth',lw);

            set(h,'color',color_nc(j,:))
            set(h(2),'LineStyle','none','Marker','s','MarkerSize',8,'markerfacecolor',color_nc(j,:))
            set(gca,'xscale','log')

            ch = get(gca,'children');
            ch = ch([4,1,2:3,5:end]);
            set(gca,'children',ch);

            drawnow
        end

    end
    set(p.h_ax(1:4),'ylim',[0.5,1])

    for i = sub_ind
        switch i
            case {1,3,4}
                dy = 0.5;
            case 2
                dy = 0.2;
        end
        Ymax = mean(rt0(i,1,:))*1.2;
        set(p.h_ax(4+i),'ylim',[0,Ymax])
        set(p.h_ax(4+i),'YTick',0:dy:Ymax);
    end
    set(p.h_ax,'xlim',[zero_coh,0.6],'xtick',[0.025,0.1,0.4],'xticklabel',100*[0.025,0.1,0.4])
    set(gcf,'OuterPosition',[100 100 1350 1050])
    p.format('FontSize',15);

    p.current_ax(1);
    ylabel('p correct')
    p.current_ax(5);
    xlabel('Coherence (%)')
    ylabel('Response time (s)')
    
    
    %% Get stats
    % This part can be used to compare mean RT between two phases 
    if 0
        for i = sub_ind
            fprintf('Subject %s Phase %s: %.4f Phase %s: %.4f\n',vSubject{i},vStage{1},mean(RT{i,1}),vStage{2},mean(RT{i,2}))
            [H(i),pval(i)] = ttest2(RT{i,1},RT{i,2});
        end

        for i = sub_ind
            for j=1:length(vStage)
                XX_lo_coh(i,j) = XX{i,j}(1);
                XX_hi_coh(i,j) = XX{i,j}(end);
            end
            rt_frac_lo_coh(i) = XX_lo_coh(i,2)/XX_lo_coh(i,1);
            rt_frac_hi_coh(i) = XX_hi_coh(i,2)/XX_hi_coh(i,1);
        end
        %%
        fprintf('RT reduction from Phase I to II\n')
        fprintf('coh =  0.0%%: %.1f %s %.1f %% RT reduction\n',(1-mean(rt_frac_lo_coh))*100,char(177),std(1-rt_frac_lo_coh)*100);
        fprintf('coh =  0.0%%: %.2f %s %.2f s to %.2f %s %.2f s\n',mean(XX_lo_coh(:,1)),char(177),std(XX_lo_coh(:,1)),mean(XX_lo_coh(:,2)),char(177),std(XX_lo_coh(:,2)));

        fprintf('coh = 51.2%%: %.1f %s %.1f %% RT reduction\n',(1-mean(rt_frac_hi_coh))*100,char(177),std(1-rt_frac_hi_coh)*100);
        fprintf('coh = 51.2%%: %.2f %s %.2f s to %.2f %s %.2f s\n',mean(XX_hi_coh(:,1)),char(177),std(XX_hi_coh(:,1)),mean(XX_hi_coh(:,2)),char(177),std(XX_hi_coh(:,2)));
    end
end


