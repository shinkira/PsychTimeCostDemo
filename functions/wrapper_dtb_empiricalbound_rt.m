function [err,P] = wrapper_dtb_empiricalbound_rt(theta,rt,coh,choice,c,pars,plot_flag,varargin)

set(0,'defaultaxesfontsize',18);
set(0,'defaulttextfontsize',18);
set(0,'defaultaxesfontweight','bold');
set(0,'defaulttextfontweight','bold');
if 1
    set(0,'defaultaxesfontsize',18);
    set(0,'defaulttextfontsize',12);
end
set(0,'defaultaxestickdir','out');
set(0,'defaultaxesbox','on');
set(0,'defaultFigureColor','w');

dt_flag = pars.dt_flag;
error_flag = pars.error_flag;

if isfield(pars,'bandwidth')
    bandwidth = pars.bandwidth;
else
    bandwidth = 0.2;
end

%%
kappa  = theta(1);
ndt_m  = theta(2);
ndt_s  = theta(3);
coh0   = theta(4);
y0a    = theta(5);

%%
dt = 0.0005;
% t  = 0:dt:6;
t  = 0:dt:10;
dy = 0.004;
y  = dy*(-(2^10):2^10)';
y0 = zeros(size(y));
y0(findclose(y,0))=1;
y0 = y0/sum(y0);

%%
% dt_method = 1; %works quite well
% dt_method = 2; %works quite well
switch dt_flag
    case 1 %ignore ndt_s for dect calculation
        dect = rt - ndt_m;
        inds = dect>0;
        q = ksdensity(dect(inds),t,'kernel','epanechnikov','support','positive','bandwidth',bandwidth,'support','positive');
        
    case 2 %deconvolve rt
        dectdist = decTimeDist(t,rt,ndt_m,ndt_s);
        q = ksdensity(t(t>0),t,'weight',dectdist(t>0),'kernel','epanechnikov','bandwidth',bandwidth,'support','positive');
        
    case 3 %add the negative mass to 0
        dect = rt - ndt_m;
        dect(dect<=0) = 0.001;
        q = ksdensity(dect,t,'kernel','epanechnikov','support','positive','bandwidth',bandwidth,'support','positive');

    case 4 %histogram
        dect = rt - ndt_m;
        edges = [0:0.025:7];
        qaux = histc(dect,edges);
        centers = (edges(1:end-1)+edges(2:end))/2;
        q = interp1(centers,qaux(1:end-1),t,'nearest');
        q(t<nanmin(centers))=0;
        q = q/nansum(q)/dt;
        
    case 5 %sskernel
        dect = rt - ndt_m;
        inds = dect>0;
        q = sskernel(dect(inds),t);

end

F = @dtb_fp_cc_searchbnd;
prior = Rtable(coh)/sum(Rtable(coh));

%%

drift = kappa * unique(coh-coh0);

dt = t(2) - t(1);
dy = y(2) - y(1);
ny = length(y);
nd = length(drift);

notabs_flag = false;
% calculating the sparse matrix outside of MEX function due to the lack of support
P = feval(F,drift,t,prior,q,y,y0,notabs_flag);

if error_flag==1
    %% eval error, binomial
    I = [choice==1, choice==0];
    n = sum(I); %all trials
    y = [sum(c(I(:,1))),sum(c(I(:,2)))]; %correct
    aux = P.up.p.*prior;
    ucoh = unique(coh);
    yfit(1) = sum(aux.*(ucoh>0)+0.5*aux.*(ucoh==0))/sum(aux);
    aux = P.lo.p.*prior;
    yfit(2) = sum(aux.*(ucoh<0)+0.5*aux.*(ucoh==0))/sum(aux);
    yfit = yfit.*n;
    err = -2 * sum(log(binopdf(y,n,yfit./n)));
    
elseif error_flag==2
    %% likelihood
    %convolve for non-decision times
    %sanity check
    if t(1)~=0
        error('for conv to work, t(1) has to be zero');
    end
    nt = length(t);
    dt = t(2)-t(1);
    ntr = length(coh);

    if ndt_s==0
        ndt = zeros(ceil(ndt_m/dt),1);
        ndt(end) = 1;
    else
        ndt = normpdf(t,ndt_m,ndt_s)*dt;
    end
    
    upRT = conv2(1,ndt(:),P.up.pdf_t);
    loRT = conv2(1,ndt(:),P.lo.pdf_t);
    upRT = upRT(:,1:nt);
    loRT = loRT(:,1:nt);
    
    rt_step = ceil(rt/dt);
    ucoh = unique(coh);
    ncoh = length(ucoh);
    p_up = nan(ntr,1);
    p_lo = nan(ntr,1);
    for i=1:ncoh
        inds = coh == ucoh(i);
        p_up(inds) = upRT(i,rt_step(inds));
        p_lo(inds) = loRT(i,rt_step(inds));
    end
    
    %clip
    p_up(p_up<eps) = eps;
    p_lo(p_lo<eps) = eps;
    pPred = p_up.*(choice==1) + p_lo.*(choice==0);
    % logl = -sum(log(pPred));
    logl = -nanmean(log(pPred));
    
    err = logl;
    
else
    err = nan;
end

%%
%% print
fprintf('err=%.3f kappa=%.2f ndt_mu=%.2f ndt_s=%.2f coh0=%.2f y0=%.2f \n',...
    err,kappa,ndt_m,ndt_s,coh0,y0a);

%%
if plot_flag
    m = prctile(rt,99.5);
    
    figure(1);clf
    
    subplot(2,2,1);
    plot(t,P.Bup,'k');
    hold all
    plot(t,P.Blo,'k');
    xlim([0,m])
    
    subplot(2,2,2);
    [tt,xx] = curva_media(choice,coh,[],0);
    plot(tt,xx,'b.-');
    hold all
    ucoh = nanunique(coh);
    plot(ucoh,P.up.p,'r-');
    legend('data','model')
    xlabel('coh')
    ylabel('p rightward')
    
    subplot(2,2,3);
    aa = sum(bsxfun(@times,P.up.pdf_t + P.lo.pdf_t, prior));
    plot(t,dt * cumsum(q),'b-');
    hold all
    plot(t,dt * cumsum(aa/dt),'r')
    %     title('cum. DT distribution')
    %     title('inexact!!')
    xlim([0,m])
    
    subplot(2,2,4);
    rt_model = (P.up.mean_t.*P.up.p+P.lo.mean_t.*P.lo.p)./(P.up.p+P.lo.p) + ndt_m; %is not exact because it
    % doesn't take into account the curtailing of the non-decision time
    % distribution !
    curva_media(rt,coh,[],1);
    hold all
    plot(ucoh,rt_model,'-');
    legend('data','model')
    xlabel('coh')
    ylabel('rt')
    
    set(gcf,'Position',[329   900  1074   415])
    format_figure(gcf);
    
    drawnow
    
end