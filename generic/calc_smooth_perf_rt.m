function out = calc_smooth_perf_rt(P,kappa,ndt_m,coh0,pars)

% addpath('/Volumes/Data/1 - LAB/15 - code/DTB_1D');

notabs_flag = false;
%kappa = P.optim.theta(1);
%ndt_m = P.optim.theta(2);
%coh0 = P.optim.theta(4);

if pars.debug
    coh = 2.^(5:9)./1e3;
else
    coh = 2.^(-2:0.1:9.2)./1e3;
end

ucoh = [fliplr(-coh),0,coh];
% ucoh = linspace(-0.6,0.6,401);
drift = kappa*(ucoh-coh0);

%run dtb
demo_flag = 0;
Pd = dtb_fp_cn_vec(drift,P.t,P.Bup,P.Blo,P.y,P.y0,notabs_flag,demo_flag);

pright = Pd.up.p;
% For Psych & Chrono plots, just compute the mean RT without the whole PDFs
% to save time and file size.
rt = ndt_m + Pd.up.p.*Pd.up.mean_t + Pd.lo.p.*Pd.lo.mean_t;

inds = drift>0;
mean_dt_c = [Pd.lo.pdf_t(~inds,:)*Pd.t'./sum(Pd.lo.pdf_t(~inds,:),2);...
    Pd.up.pdf_t(inds,:)*Pd.t'./sum(Pd.up.pdf_t(inds,:),2)];

mean_dt_nc = [Pd.up.pdf_t(~inds,:)*Pd.t'./sum(Pd.up.pdf_t(~inds,:),2);...
    Pd.lo.pdf_t(inds,:)*Pd.t'./sum(Pd.lo.pdf_t(inds,:),2)];

rt_c = mean_dt_c + ndt_m;
rt_nc = mean_dt_nc + ndt_m;

out.ucoh = ucoh;
out.pright = pright;
out.rt = rt;
out.rt_c = rt_c;
out.rt_nc = rt_nc;
