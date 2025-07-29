function opt = do_fit(fn_wrapper,dat,tl,th,tg,pars,varargin)

plot_flag = false;
fn_fit = @(theta) (fn_wrapper(theta,dat.rt,dat.coh,dat.choice,dat.c,pars,plot_flag,varargin));

if ~pars.debug
    options = optimset('Display','final','TolFun',.001);
else
    options = optimset('Display','final','TolFun',.001,'MaxFunEvals',1);
end

[theta, fval, exitflag, output] = fminsearchbnd(@(theta) fn_fit(theta),tg,tl,th,options);

opt.theta = theta;
opt.fval = fval;
opt.exitflag = exitflag;
opt.output = output;
opt.tl = tl;
opt.th = th;
opt.tg = tg;
